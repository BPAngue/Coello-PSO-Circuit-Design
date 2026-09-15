#include "psomatrixcircuit.h"
#include "random.h"

#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace PSwarm {
std::mutex consoleMutex;
std::ofstream migrationLogFile;

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

/* Migration Message where each island has one outgoing migration buffer */
struct MigrationBuffer {
    std::vector<PSwarm::Particle> particles;
};

/* Select the best particle from an Island */
std::vector<PSwarm::Particle> selectBestMigrants(const PSwarm::Swarm& swarm, unsigned migrationRate) {
    if (swarm.tPop == 0 || migrationRate == 0) {
        return {};
    }

    // number of migrants
    unsigned count = (swarm.tPop * migrationRate) / 100;

    /* If migration is enabled but population is small guarantee at least one migrant */
    if (count == 0) {
        count = 1;
    }

    count = std::min(count, swarm.tPop);

    /* Create indices */
    std::vector<unsigned> indices(swarm.tPop);

    for (unsigned i = 0; i < swarm.tPop; ++i) {
        indices[i] = i;
    }

    /* Sort best -> worst */
    std::sort(indices.begin(), indices.end(), [&](unsigned a, unsigned b) {
        return swarm.population[a].fitness > swarm.population[b].fitness;
    });

    /* Copy the best particles */
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

    /* Sort destination particles worst -> best */
    std::vector<unsigned> indices(swarm.tPop);

    for (unsigned i = 0; i < swarm.tPop; ++i) {
        indices[i] = i;
    }

    std::sort(indices.begin(), indices.end(), [&](unsigned a, unsigned b) {
        return swarm.population[a].fitness < swarm.population[b].fitness;
    });

    const unsigned count = std::min(static_cast<unsigned>(migrants.size()), swarm.tPop);

    for (unsigned i = 0; i < count; ++i) {
        const unsigned destination = indices[i];

        /* Replace the particle */
        swarm.population[destination] = migrants[i];

        /* The migrated particle becomes the particle's personal best as well */
        swarm.bestIndividualExp[destination] = migrants[i];
    }
}

/* Replace Random Particles */
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

/* Adaptive Migration Probability */
double calculateMigrationProbability(double transmitterFitness, double receiverFitness, double temperature) {
    if (temperature <= 0.0) {
        return 0.0;
    } // ask about this

    /* Receiver is better than transmitter */
    const double delta = receiverFitness - transmitterFitness;

    if (delta <= 0.0) {
        return 1.0;
    } // ask about this

    double probability = std::exp(-delta / temperature);

    /* Numerical protection */ // ask about this
    probability = std::clamp(probability, 0.0, 1.0);

    return probability;
}

/* Choose Adaptive Strategy */
MigrationStrategy chooseStrategy(double transmitterFitness, double receiverFitness, double temperature, double& probability) {
    if (transmitterFitness > receiverFitness) {
        probability = 1.0;
        return MigrationStrategy::BestWorst;
    }

    /* Otherwise calculate probability of best-random */
    probability = calculateMigrationProbability(transmitterFitness, receiverFitness, temperature);

    // Random replacement strategy
    if (rndF() < probability) {
        return MigrationStrategy::BestRandom;
    }

    return MigrationStrategy::BestWorst;
}

/* Island thread work */
void runIsland(unsigned islandIndex, PSwarm::Swarm& swarm, unsigned migrationInterval, unsigned migrationRate, 
               bool adaptive, double initialTemperature, Barrier& migrationBarrier, 
               std::vector<MigrationBuffer>& migrationBuffers, std::vector<std::unique_ptr<PSwarm::Swarm>>& islands) {
    /* Initialize this island */
    swarm.initVariables();
    swarm.reserveMemory();

    /* Unique output filename */
    swarm.nfGen = "csvs/" + swarm.nfGen + "_island" + std::to_string(islandIndex);
    swarm.nfRun = "csvs/" + swarm.nfRun;

    const std::string suffix = "_island" + std::to_string(islandIndex);

    if (swarm.nfRun.size() >= 4) {
        swarm.nfRun.insert(swarm.nfRun.size() - 4, suffix);
    } else {
        swarm.nfRun += suffix;
    }

    /* Initial Temperature */
    double temperature = initialTemperature;

    /* Independent Runs */ // this is for the number of independent runs per island not the number of generations per island
    for (unsigned run = 0; run < swarm.nRun; ++run) {
        /* Initialize run statistics */
        initStatistics(swarm.Run);

        /* Generate initial population */
        swarm.initPopulation();

        /* Per-run generation file */
        std::string fileGen;

        swarm.runFileName(run, fileGen);

        /* Seed this island. Each thread has it's own RNG state */
        const unsigned seed = initRandom(0);

        /* Run header */
        runHeader(fileGen, seed, swarm);

        {
            std::lock_guard<std::mutex> lock(consoleMutex);

            std::printf("\n[ISLAND %u] Run %u started | seed = %u", islandIndex, run, seed);
        }

        /* Generation Loop */
        for (unsigned gen = 0; gen < swarm.nGen; ++gen) {
            /* Evaluate Population */
            initStatistics(swarm.Gen);
            swarm.evaluatePopulation(gen);

            /* Record generation information */
            swarm.runInfo(fileGen, gen);

            /* Original PSO Update Coello */
            swarm.PSOAlgorithm(gen);

            /* Mutation */
            swarm.mutation();

            /* Check migration interval */
            const bool migrationTime = ((gen + 1) % migrationInterval == 0);

            if (migrationTime) {
                /* Re-evaluate new population produced by PSO + mutation */
                initStatistics(swarm.Gen);
                swarm.evaluatePopulation(gen);
                
                /* Step 1: Select emigrants */
                migrationBuffers[islandIndex].particles = selectBestMigrants(swarm, migrationRate);

                /* Everyone must finish creating their migration snapshot before anyone receives anything */
                migrationBarrier.wait();

                /* Step 2: Determine predecessor using Ring Topology */
                const unsigned numIslands = static_cast<unsigned>(islands.size());

                const unsigned transmitter = (islandIndex + numIslands - 1) % numIslands;

                /* Gen.best represents the island's best solution after the PSO update and mutation, immediately before migration. */
                const double transmitterFitness = islands[transmitter]->Gen.best.fitness;

                const double receiverFitness = swarm.Gen.best.fitness;

                /* Default strategy */
                MigrationStrategy strategy = MigrationStrategy::BestWorst;
                double probability = 0.0;

                /* check if adaptive */
                if (adaptive) {
                    strategy = chooseStrategy(transmitterFitness, receiverFitness, temperature, probability);
                }

                /* Get the sender's snapshot */
                const auto& migrants = migrationBuffers[transmitter].particles;

                /* Each thread modifies only its own swarm */
                if (strategy == MigrationStrategy::BestRandom) {
                    replaceRandom(swarm, migrants);
                } else {
                    replaceWorst(swarm, migrants);
                }

                /* Log migration */
                {
                    std::lock_guard<std::mutex> lock(consoleMutex);

                    if (adaptive) {
                        std::printf("\n[GEN %u] [MIGRATION] Island %u <- Island %u" 
                                " | Tx = %.3f | Rx = %.3f"  
                                " | pr = %.6f | %s | migrants = %zu", 
                                gen, islandIndex, transmitter, 
                                transmitterFitness, receiverFitness, 
                                probability, 
                                strategy == MigrationStrategy::BestRandom ? "best-random" : "best-worst", 
                                migrants.size());

                        if (migrationLogFile.is_open()) {
                            migrationLogFile 
                                << "[GEN " << gen << "] [MIGRATION] Island " << islandIndex
                                << " <- Island " << transmitter
                                << " | Tx = " << transmitterFitness
                                << " | Rx = " << receiverFitness
                                << " | pr = " << probability
                                << " | "
                                << (strategy == MigrationStrategy::BestRandom
                                        ? "best-random"
                                        : "best-worst")
                                << " | migrants = " << migrants.size()
                                << '\n';

                            migrationLogFile.flush();
                        }
                    } else {
                        std::printf("\n[GEN %u] [MIGRATION] Island %u <- Island %u"
                                " | best-worst | migrants = %zu",
                                gen,
                                islandIndex,
                                transmitter,
                                migrants.size());

                        if (migrationLogFile.is_open()) {
                            migrationLogFile 
                                << "[GEN " << gen << "] [MIGRATION] Island " << islandIndex
                                << " <- Island " << transmitter
                                << " | best-worst"
                                << " | migrants = " << migrants.size()
                                << '\n';

                            migrationLogFile.flush();
                        }
                    }
                }

                /* Step 3: Wait until every island has completed replacement */
                migrationBarrier.wait();

                /* Step 4: Cool Temperature */
                if (adaptive) {
                    const double alpha = rndF();
                    temperature *= alpha;

                    /* Avoid numerical zero */
                    if (temperature < 1.0e-12) {
                        temperature = 1.0e-12;
                    }
                }

                /* Step 5: Synchronize again */
                migrationBarrier.wait();
            }
        }

        /* End Run */
        runFooter(fileGen, swarm);

        {
            std::lock_guard<std::mutex> lock(consoleMutex);

            std::printf("\n[ISLAND %u] Run %u finished | Best fitness = %.3f", islandIndex, run, swarm.Run.best.fitness);
        }
    }

    /* Free memory */
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
    std::filesystem::create_directories("csvs/migration_logs");

    const std::string inputFile = argv[1];

    PSwarm::migrationLogFile.open("csvs/migration_logs/migration_log.txt", std::ios::out | std::ios::trunc);

    if (!PSwarm::migrationLogFile.is_open()) {
        std::printf("Error: Unable to open migration log file for writing.\n");
        return 1;
    }

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
        "Input file             : %s\n"
        "Number of islands      : %u\n"
        "Migration interval mu  : %u\n"
        "Migration rate chi     : %u%%\n"
        "Migration topology     : Unidirectional Ring\n"
        "Migration strategy     : %s\n"
        "Initial temperature    : %.6f\n"
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

    // Hello test
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

    /* Shared migration buffers */
    std::vector<PSwarm::MigrationBuffer> migrationBuffers(numIslands);

    /* Shared barrier */
    PSwarm::Barrier migrationBarrier(numIslands);

    /* Launch n Islands */
    std::vector<std::thread> workers;
    workers.reserve(numIslands);

    for (unsigned i = 0; i < numIslands; ++i) {
        workers.emplace_back(PSwarm::runIsland, i, std::ref(*islands[i]), migrationInterval, migrationRate, adaptive, initialTemperature, std::ref(migrationBarrier), std::ref(migrationBuffers), std::ref(islands));
    }

    /* Wait for all sub-swarms to finish */
    for (auto& worker : workers) {
        worker.join();
    }

    std::printf("\nAll %u islands finished.\n", numIslands);

    PSwarm::migrationLogFile.close();

    return 0;
}