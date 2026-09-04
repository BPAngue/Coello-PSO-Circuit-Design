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
    unsigned tNeighSwarms; // ring neighborhood size for sub-swarms

    /* Each element contains the Gen.best of one sub-swarm */
    std::vector<Particle> swarmBests; // one Gen.best per sub-swarm
    std::vector<Particle> swarmNeighborhoodBests; // one neighborhood-best per sub-swarm

    /* Synchronization point */
    Barrier barrier;

    /* Protects access to swarmBests/globalBest while they
       are being written/read */
    std::mutex bestMutex;

    ParallelState(unsigned count, unsigned neighSwarmsSize) : numSwarms(count), tNeighSwarms(neighSwarmsSize), swarmBests(count), swarmNeighborhoodBests(count), barrier(count) {}
};

/* Find the global best solution among all swarmBests */
void calculateSwarmBest(unsigned instanceIndex, Swarm& swarm, ParallelState& state) {
    {
        std::lock_guard<std::mutex> lock(state.bestMutex);
        state.swarmBests[instanceIndex] = swarm.Gen.best;
    }

    state.barrier.wait();

    if (instanceIndex == 0) {
        std::lock_guard<std::mutex> lock(state.bestMutex);

        for (unsigned i = 0; i < state.numSwarms; ++i) {
            int p = static_cast<int>(i);
            int p0 = p, p1 = p;

            for (unsigned j = 0; j < state.tNeighSwarms; j += 2) {
                p0 = (p0 + 1 < static_cast<int>(state.numSwarms)) ? p0 + 1 : 0;
                p1 = (p1 - 1 >= 0) ? p1 - 1 : static_cast<int>(state.numSwarms) - 1;

                if (state.swarmBests[static_cast<unsigned>(p0)].fitness >
                    state.swarmBests[static_cast<unsigned>(p)].fitness) p = p0;
                if (state.swarmBests[static_cast<unsigned>(p1)].fitness >
                    state.swarmBests[static_cast<unsigned>(p)].fitness) p = p1;
            }

            state.swarmNeighborhoodBests[i] = state.swarmBests[static_cast<unsigned>(p)];
        }

        // logging only
        const auto bestIt = std::max_element(
            state.swarmBests.begin(), state.swarmBests.end(),
            [](const Particle& a, const Particle& b) { return a.fitness < b.fitness; });
        std::printf("\n[SWARM] Best fitness this sync = %.2f\n", bestIt->fitness);
    }

    state.barrier.wait();
}

/* Apply globalBest to every sub-swarm */
void applyGlobalBest(unsigned instanceIndex, Swarm& swarm, ParallelState& state)
{
    const double phi3 = 1.0;
    const Particle& nBest = state.swarmNeighborhoodBests[instanceIndex];

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
                    phi3p * (static_cast<double>(nBest.chromX[d])
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
                        ? nBest.chromX[d]
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
            applyGlobalBest(instanceIndex, swarm, state);

            /* Mutation after the Global-best phase */
            swarm.mutation();

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
        std::printf("Usage: %s <input_file> [num_swarms] [M] [N] [neighborhood]\n", argv[0]);
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

    unsigned neighborhood = (argc >= 6) ? static_cast<unsigned>(std::stoul(argv[5])) : 2;

    std::printf(
        "\n====================================================\n"
        " Parallel Swarms Oriented PSO\n"
        "====================================================\n"
        "Input file       : %s\n"
        "Sub-swarms (K)   : %u\n"
        "Outer iterations : %u\n"
        "Local iterations : %u\n"
        "Neighborhood     : %u\n"
        "====================================================\n\n",
        inputFile.c_str(),
        numSwarms,
        M,
        N,
        neighborhood  
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
    PSwarm::ParallelState state(numSwarms, neighborhood);

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