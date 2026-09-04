// #include "psomatrixcircuit.h"
// #include "random.h"

// #include <string>
// #include <thread>
// #include <vector>
// #include <mutex>
// #include <atomic>
// #include <filesystem>

// namespace PSwarm {

// std::mutex consoleMutex;

// void runInstance(unsigned instanceIndex, const std::string& inputFile) {
    
//     PSwarm::Swarm swarm;

//     if (!swarm.loadParameters(inputFile)) {
//         std::lock_guard<std::mutex> lock(consoleMutex);
//         std::printf("[instance %u] failed to load parameters, aborting.\n", instanceIndex);
//         return;
//     }

//     /* Make unique output files for each instance */
//     swarm.nfGen = "csvs/" + swarm.nfGen;
//     swarm.nfGen += "_inst" + std::to_string(instanceIndex) + "_";

//     swarm.nfRun = "csvs/" + swarm.nfRun;
//     swarm.nfRun.insert(swarm.nfRun.size() - 4, "_inst" + std::to_string(instanceIndex));

//     swarm.initVariables();
//     swarm.reserveMemory();

//     globalHeader(swarm.nfRun, swarm);

//     for (unsigned run = 0; run < swarm.nRun; ++run) {
//         {
//             std::lock_guard<std::mutex> lock(consoleMutex);
//             std::printf("[instance %u] Run %02u Started\n", instanceIndex, run);
//         }

//         initStatistics(swarm.Run);
//         swarm.pSwarm(run);
//         runStatistics(swarm.nfRun, run, swarm);

//         {
//             std::lock_guard<std::mutex> lock(consoleMutex);
//             std::printf("[instance %u] Run %02u Finished\n", instanceIndex, run);
//         }
//     }

//     swarm.freeMemory();

//     std::lock_guard<std::mutex> lock(consoleMutex);
//     std::printf("[instance %u] Complete.\n", instanceIndex);
// } 

// } // namespace PSwarm


// int main(int argc, char *argv[]) {

//     if (argc < 2) {
//         std::printf("Usage: %s <input_file> [num_instances]\n", argv[0]);
//         return 0;
//     }

//     std::filesystem::create_directories("csvs");

//     const std::string inputFile = argv[1];

//     // check if arguments contain number of threads, if not use the device's max number of threads
//     unsigned numInstances = (argc >= 3) ? static_cast<unsigned>(std::stoul(argv[2])) : std::thread::hardware_concurrency();

//     // make sure that there is at least 1 thread working.
//     if (numInstances == 0) {
//         numInstances = 1;
//     }

//     std::printf("Launching %u independent solver instance(s) on '%s'\n", numInstances, inputFile.c_str());

//     std::vector<std::thread> workers;
//     workers.reserve(numInstances);

//     for (unsigned i = 0; i < numInstances; ++i) {
//         workers.emplace_back(PSwarm::runInstance, i, inputFile);
//     }

//     for (auto& t : workers) {
//         t.join();
//     }

//     std::printf("All %u instance(s) finished.\n", numInstances);

//     return 0;
// }

#include "psomatrixcircuit.h"
#include "random.h"

#include <string>
#include <thread>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <filesystem>
#include <memory>
#include <algorithm>

namespace PSwarm {

std::mutex consoleMutex;

class Barrier {
private:
    std::mutex mtx;
    std::condition_variable cv;

    unsigned threadCount;
    unsigned arrived = 0;
    unsigned generation = 0;

public:
    explicit Barrier(unsigned count) : threadCount(count) {}

    void wait()
    {
        std::unique_lock<std::mutex> lock(mtx);

        const unsigned currentGeneration = generation;

        ++arrived;

        if (arrived == threadCount) {
            arrived = 0;
            ++generation;
            cv.notify_all();
        } else {
            cv.wait(lock, [&]() {
                return generation != currentGeneration;
            });
        }
    }
};

/* Shared State between all sub-swarms */
struct ParallelState {
    unsigned numSwarms;

    /* Each element contains the Gen.best of one sub-swarm */
    std::vector<Particle> swarmBests;

    /* Best Gen.best among all sub-swarms */
    Particle globalBest;

    /* Synchronization point */
    Barrier barrier;

    /* Protects access to swarmBests/globalBest while they
       are being written/read */
    std::mutex bestMutex;

    ParallelState(unsigned count) : numSwarms(count), swarmBests(count), barrier(count) {}
};

/* Find the global best solution among all swarmBests */
void calculateSwarmBest(unsigned instanceIndex, Swarm& swarm, ParallelState& state) {
    /* Store this sub-swarm's best solution */
    {
        std::lock_guard<std::mutex> lock(state.bestMutex);

        state.swarmBests[instanceIndex] = swarm.Gen.best;
    }

    /* Wait until EVERY sub-swarm has written its Gen.best */
    state.barrier.wait();

    /* Only one thread needs to calculate the best
       We use instance 0 as the coordinator */
    if (instanceIndex == 0) {
        std::lock_guard<std::mutex> lock(state.bestMutex);
        state.globalBest = state.swarmBests[0];

        for (unsigned i = 1; i < state.numSwarms; ++i) {
            if (state.swarmBests[i].fitness > state.globalBest.fitness) {
                state.globalBest = state.swarmBests[i];
            }
        }

        std::printf("\n[SWARM] Global best fitness = %.2f\n", state.globalBest.fitness);
    }

    /* Wait until instance 0 has finished calculating globalBest */
    state.barrier.wait();
}

/* Apply globalBest to every sub-swarm */
void applyGlobalBest(Swarm& swarm, ParallelState& state)
{
    const double phi3 = 1.0;

    for (unsigned i = 0; i < swarm.tPop; ++i) {
        for (unsigned d = 0; d < swarm.nAllele; ++d) {
            const double phi1p = rndF() * swarm.phi1;
            const double phi2p = rndF() * swarm.phi2;
            const double phi3p = rndF() * phi3;

            /* Single evolutionary phase: 
               pbest = bestIndividualExp[i]
               sbest = bestSocialExp[i]
               gbest = state.globalBest */
            swarm.population[i].vi[d] +=
                    phi1p * (static_cast<double>(swarm.bestIndividualExp[i].chromX[d])
                            - static_cast<double>(swarm.population[i].chromX[d]));
            
            swarm.population[i].vi[d] += 
                    phi2p * (static_cast<double>(swarm.bestSocialExp[i].chromX[d])
                            - static_cast<double>(swarm.population[i].chromX[d]));
                            
            swarm.population[i].vi[d] += 
                    phi3p * (static_cast<double>(state.globalBest.chromX[d])
                            - static_cast<double>(swarm.population[i].chromX[d]));

            /* Clamp velocity */
            swarm.population[i].vi[d] = std::clamp(swarm.population[i].vi[d], -swarm.vMax, swarm.vMax);

            /* Update position via sigmoid-mapped velocity */
            const double vNorm = sigmoid(swarm.population[i].vi[d]);

            switch (swarm.representation) {
                case BINARY:
                    swarm.population[i].chromX[d] = flip(vNorm);
                    break;
                case INTEGER_A:
                case INTEGER_B:
                    /* do code here for integer position update */ /* Preliminary code */
                    swarm.population[i].chromX[d] = flip(vNorm)
                        ? state.globalBest.chromX[d]
                        : flip(1.0 - vNorm)
                            ? swarm.bestIndividualExp[i].chromX[d]
                            : swarm.population[i].chromX[d];
                    break;
            }
        }
    }
}

/* Run one sub-swarm */
void runInstance(unsigned instanceIndex, Swarm& swarm, ParallelState& state, unsigned M, unsigned N)
{
    /* Unique output filenames */
    swarm.nfGen = "csvs/" + swarm.nfGen;
    swarm.nfGen += "_swarm" + std::to_string(instanceIndex) + "_";

    swarm.nfRun = "csvs/" + swarm.nfRun;

    const std::string suffix = "_swarm" + std::to_string(instanceIndex);

    if (swarm.nfRun.size() >= 4) {
        swarm.nfRun.insert(swarm.nfRun.size() - 4, suffix);
    }

    /* Allocate PSO memory */
    swarm.initVariables();
    swarm.reserveMemory();

    /* One thread owns one Swarm object. */
    globalHeader(swarm.nfRun, swarm);

    /* Independent Runs */
    for (unsigned run = 0; run < swarm.nRun; ++run) {
        {
            std::lock_guard<std::mutex> lock(consoleMutex);
            std::printf("\n[SWARM %u] Run %02u Started\n", instanceIndex, run);
        }

        /* Initialize this run */
        initStatistics(swarm.Run);
        swarm.initPopulation();
        
        /* Outer LOOP M iterations */
        for (unsigned m = 0; m < M; ++m) {
            {
                std::lock_guard<std::mutex> lock(consoleMutex);

                if (instanceIndex == 0) {
                    std::printf("\n========== Outer Iteration %u / %u ==========\n", m + 1, M);
                }
            }

            /* Local Independent PSO Evolution N iterations */
            for (unsigned n = 0; n < N; ++n) {
                const unsigned generation = m * N + n; // n = 0 is not use because evaluatePopulation(0) means initialize/update personal best
                initStatistics(swarm.Gen);
                swarm.evaluatePopulation(generation);
                
                /* Update normal run statistics */
                swarm.runInfo("", generation);

                /* Normal LOCAL PSO update */
                swarm.PSOAlgorithm(generation);
                swarm.mutation();
            }

            const unsigned globalGeneration = (m + 1) * N;

            swarm.evaluatePopulation(globalGeneration);

            /* SYNCHRONIZATION POINT */
            calculateSwarmBest(instanceIndex, swarm, state);

            /* GLOBAL-BEST PHASE */
            applyGlobalBest(swarm, state);

            /* Mutation after the Global-best phase */
            // swarm.mutation();

            /* Wait for every swarm to finish */
            state.barrier.wait();
        }

        /* END OF RUN */
        { 
            std::lock_guard<std::mutex> lock(consoleMutex);
            std::printf("\n[Swarm %u] Run %02u Finished\n", instanceIndex, run);
        }
    }

    swarm.freeMemory();

    {
        std::lock_guard<std::mutex> lock(consoleMutex);

        std::printf("[Swarm %u] Complete.\n", instanceIndex);
    }
}

} // namespace PSwarm

/* MAIN ENTRY POINT */
int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::printf("Usage: %s <input_file> [num_swarms] [M] [N]\n", argv[0]);
        return 0;
    }

    std::filesystem::create_directories("csvs");

    const std::string inputFile = argv[1];

    /* K: number of parallel sub-swarms */
    unsigned numSwarms = (argc >= 3) ? static_cast<unsigned>(std::stoul(argv[2])) : std::thread::hardware_concurrency();

    if (numSwarms == 0) {
        numSwarms = 1;
    }

    /* M: number of outer parallel-swarm iterations 
       DEFAULT: 10 */
    unsigned M = (argc >= 4) ? static_cast<unsigned>(std::stoul(argv[3])) : 10;

    /* N: number of independent local PSO generations by each subswarm 
       DEFAULT: 5 */
    unsigned N = (argc >= 5) ? static_cast<unsigned>(std::stoul(argv[4])) : 5;

    std::printf(
        "\n====================================================\n"
        " Parallel Swarms Oriented PSO\n"
        "====================================================\n"
        "Input file       : %s\n"
        "Sub-swarms (K)   : %u\n"
        "Outer iterations : %u\n"
        "Local iterations : %u\n"
        "====================================================\n\n",
        inputFile.c_str(),
        numSwarms,
        M,
        N  
    );

    /* Create K independent Swarm objects. 
       Each object represents one sub-swarm. */
    std::vector<std::unique_ptr<PSwarm::Swarm>> swarms;

    swarms.reserve(numSwarms);

    for (unsigned i = 0; i < numSwarms; ++i) {
        auto swarm = std::make_unique<PSwarm::Swarm>();

        if (!swarm->loadParameters(inputFile)) {
            std::printf("Failed to load parameters for swarm %u.\n", i);
            return 1;
        }

        swarms.push_back(std::move(swarm));
    }

    /* Shared synchronizaiton / communication state */
    PSwarm::ParallelState state(numSwarms);

    /* Launch K sub-swarms */
    std::vector<std::thread> workers;
    workers.reserve(numSwarms);

    for (unsigned i = 0; i < numSwarms; ++i) {
        workers.emplace_back(PSwarm::runInstance, i, std::ref(*swarms[i]), std::ref(state), M, N);
    }

    /* Wait for all sub-swarms to finish */
    for (auto& worker : workers) {
        worker.join();
    }

    std::printf("\nAll %u sub-swarms finished.\n", numSwarms);

    return 0;
}