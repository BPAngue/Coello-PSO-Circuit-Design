/********************************************************/
/*                  CINVESTAV - IPN                     */
/*        Department of Electrical Engineering          */
/*                 Computing Section                    */
/*                                                      */
/*               Evolutionary Computation               */
/*                                                      */
/* File: psomatrixcircuit_mt.cpp                        */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Multithreaded entry point. Spawns N      */
/* worker threads, each owning a completely independent  */
/* PSwarm::Swarm instance (own population, own decoder   */
/* buffers, own thread-local RNG stream, own output      */
/* files). This is "N separate solver instances running  */
/* concurrently" -- NOT population-level parallelism      */
/* inside a single Swarm.                                */
/*                                                       */
/* Usage:                                                */
/*   psomatrixcircuit_mt <input_file> [num_instances]    */
/*                                                       */
/* If num_instances is omitted, it defaults to           */
/* std::thread::hardware_concurrency().                  */
/*                                                       */
/* IMPORTANT PRECONDITIONS (see accompanying notes):     */
/*  - random.h/.cpp must use a thread_local RNG engine   */
/*    (patched -- std::rand()/srand() are process-global */
/*    and are a data race across threads).                */
/*  - Each instance gets its own output filenames         */
/*    (nfRun/nfGen suffixed with the instance index) so   */
/*    concurrent file writes in statistics.cpp never      */
/*    collide.                                            */
/********************************************************/

#include "psomatrixcircuit.h"
#include "random.h"

#include <cstdio>
#include <string>
#include <thread>
#include <vector>
#include <mutex>
#include <atomic>
#include <random>

namespace {

/* Guards only the human-readable progress prints below so
   lines from different threads don't interleave mid-sentence.
   It does NOT guard any library code -- statistics.cpp writes
   to per-instance files, so no locking is needed there. */
std::mutex consoleMutex;

/* Runs one completely independent solver instance on the
   calling thread. Mirrors the single-threaded main() in
   psomatrixcircuit.cpp, but:
     - suffixes output filenames with the instance index
     - seeds this thread's RNG from a distinct value derived
       from the instance index, so instances don't correlate
       even if launched in the same wall-clock instant */
void runInstance(const std::string& inputFile, unsigned instanceIndex, unsigned baseSeed)
{
    PSwarm::Swarm swarm;

    if (!swarm.loadParameters(inputFile)) {
        std::lock_guard<std::mutex> lock(consoleMutex);
        std::printf("[instance %u] failed to load parameters, aborting.\n", instanceIndex);
        return;
    }

    /* Make this instance's output files unique so concurrent
       writers never touch the same path. */
    swarm.nfGen += "_inst" + std::to_string(instanceIndex) + "_";
    swarm.nfRun.insert(swarm.nfRun.size() - 4 /* before ".csv" */,
                        "_inst" + std::to_string(instanceIndex));

    swarm.initVariables();
    swarm.reserveMemory();

    /* Give this instance a distinct, deterministic base seed.
       pSwarm() combines this with the run index each time it's
       called, so no two (instance, run) pairs anywhere ever
       produce the same seed, regardless of timing. */
    swarm.instanceSeed = baseSeed + instanceIndex * 1000u; // spacing avoids overlap even across many runs

    globalHeader(swarm.nfRun, swarm);

    for (unsigned run = 0; run < swarm.nRun; ++run) {
        {
            std::lock_guard<std::mutex> lock(consoleMutex);
            std::printf("[instance %u] Run %02u Started\n", instanceIndex, run);
        }

        initStatistics(swarm.Run);
        swarm.pSwarm(run);
        runStatistics(swarm.nfRun, run, swarm);

        {
            std::lock_guard<std::mutex> lock(consoleMutex);
            std::printf("[instance %u] Run %02u Finished\n", instanceIndex, run);
        }
    }

    swarm.freeMemory();

    std::lock_guard<std::mutex> lock(consoleMutex);
    std::printf("[instance %u] Complete.\n", instanceIndex);
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::printf("Usage: %s <input_file> [num_instances]\n", argv[0]);
        return 0;
    }

    const std::string inputFile = argv[1];

    unsigned numInstances = (argc >= 3)
        ? static_cast<unsigned>(std::stoul(argv[2]))
        : std::thread::hardware_concurrency();

    if (numInstances == 0) numInstances = 1;

    std::printf("Launching %u independent solver instance(s) on '%s'\n",
                numInstances, inputFile.c_str());

    /* A base seed shared across instances, offset per instance
       inside runInstance(). Using random_device once here (main
       thread) rather than time(nullptr) avoids every run using
       the same seed if the program is launched repeatedly within
       the same second. */
    const unsigned baseSeed = std::random_device{}();

    std::vector<std::thread> workers;
    workers.reserve(numInstances);

    for (unsigned i = 0; i < numInstances; ++i)
        workers.emplace_back(runInstance, inputFile, i, baseSeed);

    for (auto& t : workers)
        t.join();

    std::printf("All %u instance(s) finished.\n", numInstances);
    return 1;
}
