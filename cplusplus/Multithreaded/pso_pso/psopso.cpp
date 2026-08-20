/********************************************************/
/*                  CINVESTAV - IPN                     */
/*        Department of Electrical Engineering          */
/*                 Computing Section                    */
/*                                                      */
/*               Evolutionary Computation               */
/*                                                      */
/*                    PSO-PSO orchestrator               */
/*                     April 2026                       */
/*                                                      */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* File: psopso.cpp                                     */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* See psopso.h for the algorithm this implements.       */
/*                                                       */
/* Threading model: for every cycle, K worker threads are */
/* spawned for the multi-evolutionary phase and joined     */
/* before the global-best is computed (the join is the      */
/* synchronisation barrier -- no swarm reads another          */
/* swarm's state while any thread is running). K threads are  */
/* then spawned again for the single-evolutionary phase and    */
/* joined the same way. gBest is only ever written by the main */
/* thread, between phases, while no worker thread is running.   */
/* Each worker seeds its own thread-local RNG stream (see        */
/* random.h) as the first thing it does, so streams never        */
/* correlate across swarms/cycles/phases.                         */
/********************************************************/

#include "psopso.h"
#include "psomatrixcircuit.h"
#include "random.h"

#include <cstdarg>
#include <cstdio>
#include <limits>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <vector>

namespace PSwarm {

namespace {

std::mutex consoleMutex;

void logf(const char* fmt, ...)
{
    std::va_list args;
    va_start(args, fmt);
    std::lock_guard<std::mutex> lock(consoleMutex);
    std::vprintf(fmt, args);
    va_end(args);
}

void appendCycleLog(const std::string& path, unsigned cycle,
                     const char* phase, double fitness)
{
    std::lock_guard<std::mutex> lock(consoleMutex);
    std::FILE* pf = std::fopen(path.c_str(), "a");
    if (!pf) return;
    std::fprintf(pf, "%u,%s,%.8f\n", cycle, phase, fitness);
    std::fclose(pf);
}

} // namespace

void runPSOPSO(const PSOPSOConfig& cfg)
{
    const unsigned K = cfg.numSwarms ? cfg.numSwarms : 1;

    std::vector<Swarm>   swarms(K);
    std::vector<unsigned> genCounter(K, 0);
    std::vector<char>     everInitialised(K, 0); /* vector<char>, NOT vector<bool> --
                                                     vector<bool> bit-packs its storage,
                                                     so concurrent writes to different
                                                     indices from different threads would
                                                     be a genuine data race. */

    /* ---- Load & prepare every sub-swarm (sequential: file I/O
       and stdout printing here, no need to parallelise it) ---- */
    for (unsigned k = 0; k < K; ++k) {
        if (!swarms[k].loadParameters(cfg.inputFile)) {
            std::printf("[PSO-PSO] sub-swarm %u failed to load '%s', aborting.\n",
                        k, cfg.inputFile.c_str());
            return;
        }
        swarms[k].initVariables();
        swarms[k].reserveMemory();
        swarms[k].initPopulation();
        swarms[k].initSwarmBest();

        /* Unique per-swarm output filenames, same convention as
           psomatrixcircuit_mt.cpp's per-instance suffixing. */
        swarms[k].nfGen += "_swarm" + std::to_string(k) + "_";
        swarms[k].nfRun.insert(swarms[k].nfRun.size() - 4 /* before ".csv" */,
                               "_swarm" + std::to_string(k));
        globalHeader(swarms[k].nfRun, swarms[k]);
    }

    const unsigned baseSeed = cfg.baseSeed ? cfg.baseSeed : std::random_device{}();

    Particle gBest;
    gBest.fitness = -std::numeric_limits<double>::infinity();

    const std::string cycleLogPath = cfg.inputFile + "_psopso_cycles.csv";
    {
        std::FILE* pf = std::fopen(cycleLogPath.c_str(), "w");
        if (pf) {
            std::fprintf(pf, "Cycle,Phase,GlobalBestFitness\n");
            std::fclose(pf);
        }
    }

    std::printf("[PSO-PSO] K=%u sub-swarms, N=%u (multi-phase), M=%u (single-phase), "
               "%u cycle(s), phi3=%.3f\n",
               K, cfg.multiPhaseIters, cfg.singlePhaseIters, cfg.numCycles, cfg.phi3);

    for (unsigned cycle = 0; cycle < cfg.numCycles; ++cycle) {

        /* =========== Multi-evolutionary phase =========== */
        {
            std::vector<std::thread> workers;
            workers.reserve(K);

            for (unsigned k = 0; k < K; ++k) {
                workers.emplace_back([&, k]() {
                    initRandom(baseSeed + k * 7919u + cycle * 104729u + 1u);

                    std::string genFile = swarms[k].nfGen + "c" + std::to_string(cycle) + "_multi.csv";
                    runHeader(genFile, baseSeed + k, swarms[k]);

                    for (unsigned g = 0; g < cfg.multiPhaseIters; ++g) {
                        swarms[k].evaluatePopulationPSOPSO(genCounter[k]);
                        generationStats(genFile, genCounter[k], swarms[k]);

                        swarms[k].PSOAlgorithmPSOPSO(genCounter[k], !everInitialised[k],
                                                     /*singlePhase=*/false, gBest, cfg.phi3);
                        swarms[k].mutation();

                        everInitialised[k] = 1;
                        ++genCounter[k];
                    }
                });
            }
            for (auto& t : workers) t.join();
        }

        /* ---- Step 6: gbest = best sbest across all K swarms ---- */
        for (unsigned k = 0; k < K; ++k)
            if (swarms[k].sBest.fitness > gBest.fitness)
                gBest = swarms[k].sBest;

        appendCycleLog(cycleLogPath, cycle, "multi", gBest.fitness);
        logf("[PSO-PSO] cycle %u/%u  after multi-phase   gbest fitness = %.6f\n",
             cycle + 1, cfg.numCycles, gBest.fitness);

        /* =========== Single-evolutionary phase =========== */
        {
            std::vector<std::thread> workers;
            workers.reserve(K);

            for (unsigned k = 0; k < K; ++k) {
                workers.emplace_back([&, k]() {
                    initRandom(baseSeed + k * 7919u + cycle * 104729u + 50000u);

                    std::string genFile = swarms[k].nfGen + "c" + std::to_string(cycle) + "_single.csv";

                    for (unsigned g = 0; g < cfg.singlePhaseIters; ++g) {
                        swarms[k].evaluatePopulationPSOPSO(genCounter[k]);
                        generationStats(genFile, genCounter[k], swarms[k]);

                        swarms[k].PSOAlgorithmPSOPSO(genCounter[k], /*isFirstGeneration=*/false,
                                                     /*singlePhase=*/true, gBest, cfg.phi3);
                        swarms[k].mutation();

                        ++genCounter[k];
                    }
                });
            }
            for (auto& t : workers) t.join();
        }

        /* Sub-swarms kept improving during the single-phase too;
           refresh gbest before it's used as next cycle's guide. */
        for (unsigned k = 0; k < K; ++k)
            if (swarms[k].sBest.fitness > gBest.fitness)
                gBest = swarms[k].sBest;

        appendCycleLog(cycleLogPath, cycle, "single", gBest.fitness);
        logf("[PSO-PSO] cycle %u/%u  after single-phase  gbest fitness = %.6f\n",
             cycle + 1, cfg.numCycles, gBest.fitness);
    }

    /* ---- Final report. Any swarm's decoder can render gBest's
       Boolean expression, since every sub-swarm was loaded from
       the same input file and therefore shares identical circuit
       dimensions / gate set. ---- */
    std::printf("\n[PSO-PSO] Finished after %u cycle(s). Best fitness found: %.8f\n",
               cfg.numCycles, gBest.fitness);
    std::printf("Violations: %d   Gates used: %u\n",
               static_cast<int>(swarms[0].circuit.numTotalOutputs)
                   - static_cast<int>(gBest.numEqual),
               gBest.numGates);
    std::printf("Boolean Expression:\n");
    for (unsigned i = 0; i < swarms[0].circuit.numOutputs; ++i) {
        std::string expr = swarms[0].decoder.expression(gBest.chromX, i, swarms[0], swarms[0].circuit);
        std::printf("%s\n", expr.c_str());
    }

    for (unsigned k = 0; k < K; ++k)
        swarms[k].freeMemory();
}

} // namespace PSwarm
