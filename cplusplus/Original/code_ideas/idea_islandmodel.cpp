#include "psomatrixcircuit.h"
#include "random.h"

#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace PSwarm {

std::mutex consoleMutex;
std::mutex statisticsMutex;

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

/* Migration Strategy Enum */
enum class MigrationStrategy {
    BestWorst,
    BestRandom
};

/* Contains the migrants (particles) that we want to emigrate into another Island */
struct MigrationMessage {
    unsigned sourceIsland;
    std::vector<PSwarm::Particle> migrants;

    MigrationMessage() : sourceIsland(0) {}
};

/* Select the best particle from an Island */
std::vector<PSwarm::Particle> selectBestMigrants(const PSwarm::Swarm& swarm, unsigned count) {
    count = std::min(count, swarm.tPop);

    std::vector<unsigned> indices(swarm.tPop);

    for (unsigned i = 0; i < swarm.tPop; ++i) {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](unsigned a, unsigned b) {
        return swarm.population[a].fitness > swarm.population[b].fitness;
    });

    std::vector<PSwarm::Particle> migrants;

    migrants.reserve(count);

    for (unsigned i = 0; i < count; ++i) {
        migrants.push_back(swarm.population[indices[i]]);
    }

    return migrants;
}

/* Replace Worst Particles */
void replaceWorst(PSwarm::Swarm& swarm, const std::vector<PSwarm::Particle>& migrants) {
    if (migrants.empty()) {
        return;
    }

    std::vector<unsigned> indices(swarm.tPop);

    for (unsigned i = 0; i < swarm.tPop; ++i) {
        indices[i] = i;
    }

    /* Sort from worst to best */
    std::sort(indices.begin(), indices.end(), [&](unsigned a, unsigned b) {
        return swarm.population[a].fitness < swarm.population[b].fitness;
    });

    const unsigned count = std::min(static_cast<unsigned>(migrants.size()), swarm.tPop);

    for (unsigned i = 0; i < count; ++i) {
        const unsigned destination = indices[i];

        /* Insert the complete migrant */
        swarm.population[destination] = migrants[i];

        /* The migrated particle becomes the particle's personal best as well */
        swarm.bestIndividualExp[destination] = migrants[i];
    }
}

void replaceRandom(PSwarm::Swarm& swarm, const std::vector<PSwarm::Particle>& migrants) {
    if (migrants.empty()) {
        return;
    }

    std::vector<unsigned> indices(swarm.tPop);

    for (unsigned i = 0; i < swarm.tPop; ++i) {
        indices[i] = i;
    }

    /* Shuffle indices to introduce randomness */
    shuffle(indices);

    const unsigned count = std::min(static_cast<unsigned>(migrants.size()), swarm.tPop);
    
    for (unsigned i = 0; i < count; ++i) {
        const unsigned destination = indices[i];

        swarm.population[destination] = migrants[i];

        /* The migrated particle becomes the particle's personal best as well */
        swarm.bestIndividualExp[destination] = migrants[i];
    }
}

double adaptiveProbability(double transmitterFitness, double receiverFitness, double temperature) {
    if (temperature <= 0.0) {
        return 0.0;
    }

    /* Difference between the receiver's better fitness and the transmitter's worse fitness */
    const double delta = receiverFitness - transmitterFitness;

    if (delta <= 0.0) {
        return 1.0;
    }

    double probability = std::exp(-delta / temperature);

    return probability;
}

/* Determine migration strategy */
MigrationStrategy chooseAdaptiveStrategy(double transmitterFitness, double receiverFitness, double temperature, double& probability) {
    if (transmitterFitness > receiverFitness) {
        probability = 1.0;
        return MigrationStrategy::BestWorst;
    }

    probability = adaptiveProbability(transmitterFitness, receiverFitness, temperature);

    if (rndF() < probability) {
        return MigrationStrategy::BestRandom;
    }

    return MigrationStrategy::BestWorst;
}

/* Run one Migration epoch */
void performRingMigration(std::vector<std::unique_ptr<PSwarm::Swarm>>& islands, unsigned migrationRate, double& temperature, bool adaptive) {
    const unsigned numIslands = static_cast<unsigned>(islands.size());

    if (numIslands < 2) {
        return;
    }

    /* Snapshot all emigrants */
    std::vector<std::vector<PSwarm::Particle>> emigrants(numIslands);

    for (unsigned i = 0; i < numIslands; ++i) {
        unsigned count = islands[i]->tPop * migrationRate / 100;

        /* Ensure at least one migrant if migration rate > 0 */
        if (migrationRate > 0 && count == 0) {
            count = 1;
        }

        count = std::min(count, islands[i]->tPop);

        emigrants[i] = selectBestMigrants(*islands[i], count);
    }

    /* Apply ring migration */
    for (unsigned receiver = 0; receiver < numIslands; ++receiver) {
        /* Previous island send to this island */
        const unsigned transmitter = (receiver + numIslands - 1) % numIslands;

        const double transmitterFitness = islands[transmitter]->Gen.best.fitness;
        const double receiverFitness = islands[receiver]->Gen.best.fitness;

        double probability = 0.0;
        
        MigrationStrategy strategy = MigrationStrategy::BestWorst;

        if (adaptive) {
            strategy = chooseAdaptiveStrategy(transmitterFitness, receiverFitness, temperature, probability);
        }

        /* Logging */
        {
            std::lock_guard<std::mutex> lock(consoleMutex);
            if (adaptive) {
                std::printf("\n[MIGRATION] Island %u -> Island %u" 
                            " | Tx gbest = %.3f" 
                            " | Rx gbest = %.3f" 
                            " | Prob = %.6f "
                            " | Strategy = %s", 
                            transmitter, 
                            receiver, 
                            transmitterFitness,
                            receiverFitness, 
                            probability,
                            strategy == MigrationStrategy::BestWorst ? "best-worst" : "best-random"
                ); 
            } else {
                std::printf("\n[MIGRATION] Island %u -> Island %u" 
                            " | Strategy = best-worst", 
                            transmitter, 
                            receiver
                );
            }
        }

        /* Apply replacement */
        if (strategy == MigrationStrategy::BestRandom) {
            replaceRandom(*islands[receiver], emigrants[transmitter]);
        } else {
            replaceWorst(*islands[receiver], emigrants[transmitter]);
        }
    }

    /* Cool down the temperature */
    if (adaptive) {
        const double alpha = rndF();
        temperature *= alpha;

         /*
         * Avoid temperature becoming exactly zero.
         */
        if (temperature < 1.0e-12) {
            temperature = 1.0e-12;
        }

        {
            std::lock_guard<std::mutex> lock(consoleMutex);

            std::printf("\n[ADAPTIVE] alpha = %.6f"
                        " | new Tem = %.10f",
                        alpha,
                        temperature
            );
        }
    }
}

void runInstance(unsigned islandIndex, PSwarm::Swarm& swarm, unsigned migrationInterval, 
                 unsigned migrationRate, bool adaptive, double initialTemperature, 
                 Barrier& barrier, std::vector<std::unique_ptr<PSwarm::Swarm>>& islands) 
{
    /* Give every island a different output filename. */

    swarm.nfGen = "csvs/" + swarm.nfGen + "_island" + std::to_string(islandIndex);

    swarm.nfRun = "csvs/" + swarm.nfRun;

    const std::string suffix = "_island" + std::to_string(islandIndex);

    if (swarm.nfRun.size() >= 4) {
        swarm.nfRun.insert(swarm.nfRun.size() - 4, suffix);
    } else {
        swarm.nfRun += suffix;
    }

    /* Allocate PSO memory */
    swarm.initVariables();
    swarm.reserveMemory();

    /* Each thread receives its own RNG. */
    const unsigned seed = initRandom(0);

    /* Temperature belongs to this island/thread. */
    double temperature = initialTemperature;

    /* Independent Runs */ // this is for the number of independent runs per island not the number of generations per island
    for (unsigned run = 0; run < swarm.nRun; ++run) {
        {
            std::lock_guard<std::mutex> lock(consoleMutex);
            std::printf("\n[ISLAND %u] Run %02u Started | Seed %u\n", islandIndex, run, seed);
        }

        /* Reset run statistics */
        initStatistics(swarm.Run);

        /* Create initial population */
        swarm.initPopulation();

        /* Create the per run filename */
        std::string fileGen;
        swarm.runFileName(run, fileGen);

        /* Write per-run header */
        runHeader(fileGen, seed, swarm);

        /* Generation loop */
        for (unsigned gen = 0; gen < swarm.nGen; ++gen) {
            /* reset generation statistics */
            initStatistics(swarm.Gen);

            /* evaluate current population */
            swarm.evaluatePopulation(gen);

            /* Record generation information */
            swarm.runInfo(fileGen, gen);

            /* Synchronize before migration */
            barrier.wait();

            /* Migration */
            if ((gen + 1) % migrationInterval == 0) {
                if (islandIndex == 0) {

                    performRingMigration(
                        islands,
                        migrationRate,
                        temperature,
                        adaptive
                    );
                }
            }

            /* Wait until migration has finished */
            barrier.wait();

            swarm.PSOAlgorithm(gen);

            /* Mutation from the original algorithm */
            swarm.mutation();
        }

        /* Finish this run */
        runFooter(fileGen, swarm);

        /* Global statistics file is shared by all islands. */
        {
            std::lock_guard<std::mutex> lock(statisticsMutex);

            runStatistics(swarm.nfRun, run, swarm);
        }

        {
            std::lock_guard<std::mutex> lock(consoleMutex);

            std::printf("\n[ISLAND %u] Run %02u finished", islandIndex, run);
        }
    }

    /* Release memory owned by this island. */
    swarm.freeMemory();
}

} // namespace PSwarm

/* MAIN ENTRY POINT */
int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::printf("Usage: %s <input_file> [num_islands] [migration_interval] [migration_rate] [adaptive] [initial_temperature]\n", argv[0]);
        return 0;
    }

    std::filesystem::create_directories("csvs");
    const std::string inputFile = argv[1];

    /* n Number of islands */
    unsigned numIslands = (argc >= 3) ? static_cast<unsigned>(std::stoul(argv[2])) : 10;

    /* μ Migration Interval */
    unsigned migrationInterval = (argc >= 4) ? static_cast<unsigned>(std::stoul(argv[3])) : 100;

    /* χ Migration Rate */
    unsigned migrationRate = (argc >= 5) ? static_cast<unsigned>(std::stoul(argv[4])) : 10;

    /* Adaptive Strategy boolean */
    bool adaptive = (argc >= 6) ? (std::stoi(argv[5]) != 0) : false;

    /* Initial Temperature */
    double initialTemperature = (argc >= 7) ? std::stod(argv[6]) : 10.0;

    /* Validate parameters. */

    if (numIslands < 2) {
        std::printf("Error: num_islands must be >= 2.\n");
        return 1;
    }

    if (migrationInterval == 0) {
        std::printf("Error: migration_interval must be > 0.\n");
        return 1;
    }

    if (migrationRate == 0 || migrationRate > 100) {
        std::printf("Error: migration_rate must be between 1 and 100.\n");
        return 1;
    }

    if (initialTemperature <= 0.0) {
        std::printf("Error: initial_temperature must be > 0.\n");
        return 1;
    }


    std::printf(
        "\n"
        "====================================================\n"
        "Island Model Particle Swarm Optimization\n"
        "====================================================\n"
        "Input file            : %s\n"
        "Number of islands     : %u\n"
        "Migration interval μ  : %u\n"
        "Migration rate χ      : %u%%\n"
        "Migration topology    : Unidirectional Ring\n"
        "Migration strategy    : %s\n"
        "Initial temperature   : %.6f\n"
        "====================================================\n\n",
        inputFile.c_str(),
        numIslands,
        migrationInterval,
        migrationRate,
        adaptive
            ? "Adaptive"
            : "Best-Worst",
        initialTemperature
    );

    /* Create n Independent Swarm objects.
       Each object represents one island. */
    std::vector<std::unique_ptr<PSwarm::Swarm>> islands;

    islands.reserve(numIslands);

    for (unsigned i = 0; i < numIslands; ++i) {
        auto swarm = std::make_unique<PSwarm::Swarm>();

        if (!swarm->loadParameters(inputFile)) {
            std::printf("Failed to load parameters for island %u.\n", i);
            return 1;
        }

        islands.push_back(std::move(swarm));
    }

    /* Initialize one global CSV header */
    globalHeader("csvs/impsp_global.csv", *islands[0]);

    /* Barrier shared by all island threads */
    PSwarm::Barrier barrier(numIslands);

    /* Launch n Islands */
    std::vector<std::thread> workers;
    workers.reserve(numIslands);

    for (unsigned i = 0; i < numIslands; ++i) {
        workers.emplace_back(PSwarm::runInstance, i, std::ref(*islands[i]), migrationInterval, migrationRate, adaptive, initialTemperature, std::ref(barrier), std::ref(islands));
    }

    /* Wait for all sub-swarms to finish */
    for (auto& worker : workers) {
        worker.join();
    }

    std::printf("\nAll %u islands finished.\n", numIslands);

    return 0;
}