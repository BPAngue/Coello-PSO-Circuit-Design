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

/* Migrants */
struct Migrant {
    Particle particle;
    Particle personalBest;
    unsigned origin = 0;
    double fitnessBefore = 0.0;
};

/* Tracking information for particles currently residing on an island */
struct ParticleTracking {
    std::vector<unsigned> origin; /* island from which population[i] originated */
    std::vector<double> fitnessBefore; /* fitness of that particle immediately before it entered this island */

    void initialize(unsigned island, const Swarm& swarm) {
        const unsigned size = static_cast<unsigned>(swarm.population.size());
        origin.assign(size, island);
        fitnessBefore.resize(size);
        for (unsigned i = 0; i < size; ++i) {
            fitnessBefore[i] = swarm.population[i].fitness;
        }
    }

    void resize(unsigned size, unsigned defaultIsland) {
        origin.resize(size, defaultIsland);
        fitnessBefore.resize(size, 0.0);
    }
};

class DIMPSOController {
public:
    unsigned nIslands = 0;
    double alpha = 0.8;
    double beta = 0.01;

    /* Migration matrix M[i][j] -> the probability of moving from i to j */
    std::vector<std::vector<double>> migrationMatrix;

    /* Data matrix Dm[i][j] -> improvement of particles from i after processing on j */
    std::vector<std::vector<double>> dataMatrix;

    /* Feedback accumulated during the current migration epoch */
    std::vector<std::vector<double>> feedbackSum;

    std::vector<std::vector<unsigned>> feedbackCount;

    /* Outgoing migrants outgoing[source][destination]*/
    std::vector<std::vector<std::vector<Migrant>>> outgoing;

    DIMPSOController() = default;

    DIMPSOController(unsigned islands, double a, double b) {
        initialize(islands, a, b);
    }

    /* Initialize DIMPSO */
    void initialize(unsigned islands, double a, double b) {
        nIslands = islands;
        alpha = a;
        beta = b;
        migrationMatrix.assign(nIslands, std::vector<double>(nIslands, 0.0));
        dataMatrix.assign(nIslands, std::vector<double>(nIslands, 0.0));
        feedbackSum.assign(nIslands, std::vector<double>(nIslands, 0.0));
        feedbackCount.assign(nIslands, std::vector<unsigned>(nIslands, 0));
        outgoing.resize(nIslands);

        for (unsigned i = 0; i < nIslands; ++i) {
            outgoing[i].resize(nIslands);
        }

        /* Fully connected initialization */
        const double initialProbability = 1.0 / static_cast<double>(nIslands);

        for (unsigned i = 0; i < nIslands; ++i) {
            for (unsigned j = 0; j < nIslands; ++j) {
                migrationMatrix[i][j] = initialProbability;
            }
        }
    }

    /* Reset for a new independent PSO run */
    void reset() {
        const double initialProbability = 1.0 / static_cast<double>(nIslands);

        for (unsigned i = 0; i < nIslands; ++i) {
            for (unsigned j = 0; j < nIslands; ++j) {
                migrationMatrix[i][j] = initialProbability;
                dataMatrix[i][j] = 0.0;
                feedbackSum[i][j] = 0.0;
                feedbackCount[i][j] = 0;
                outgoing[i][j].clear();
            }
        }
    }

    /* Clear feedback produced by analyze() */
    void clearFeedback() {
        for (unsigned i = 0; i < nIslands; ++i) {
            for (unsigned j = 0; j < nIslands; ++j) {
                feedbackSum[i][j] = 0.0;
                feedbackCount[i][j] = 0;
            }
        }
    }

    void analyze(unsigned destination, const Swarm& swarm, const ParticleTracking& tracking) {
        const unsigned size = static_cast<unsigned>(swarm.population.size());

        for (unsigned p = 0; p < size; ++p) {
            const unsigned source = tracking.origin[p];

            if (source >= nIslands) {
                continue;
            }

            const double improvement = swarm.population[p].fitness - tracking.fitnessBefore[p];

            feedbackSum[source][destination] += improvement;
            feedbackCount[source][destination] += 1;
        }
    }

    /* Learn */
    void learn() {
        for (unsigned i = 0; i < nIslands; ++i) {
            for (unsigned j = 0; j < nIslands; ++j) {
                if (feedbackCount[i][j] == 0) {
                    continue;
                }

                dataMatrix[i][j] = feedbackSum[i][j] / static_cast<double>(feedbackCount[i][j]);
            }
        }
    }

    /* Generate stochastic vector Y */
    std::vector<double> generateStochasticVector() const {
        std::vector<double> Y(nIslands, 0.0);
        double total = 0.0;

        for (unsigned j = 0; j < nIslands; ++j) {
            Y[j] = rndF() + 1.0e-12;
            total += Y[j];
        }

        if (total <= 0.0) {
            const double p = 1.0 / static_cast<double>(nIslands);
            std::fill(Y.begin(), Y.end(), p);

            return Y;
        } 

        for (double& value : Y) {
            value /= total;
        }

        return Y;
    }

    void updateMigrationMatrix(unsigned island) {
        if (nIslands == 0) {
            return;
        }

        /* Find B = argmax Dm(i,j) */
        double maximum = dataMatrix[island][0];

        for (unsigned j = 1; j < nIslands; ++j) {
            if (dataMatrix[island][j] > maximum) {
                maximum = dataMatrix[island][j];
            }
        }

        /* Reward vector R */
        std::vector<unsigned> B;

        const double epsilon = 1.0e-12;

        for (unsigned j = 0; j < nIslands; ++j) {
            if (std::fabs(dataMatrix[island][j] - maximum) <= epsilon) {
                B.push_back(j);
            }
        }

        std::vector<double> R(nIslands, 0.0);

        if (!B.empty()) {
            const double reward = 1.0 / static_cast<double>(B.size());
            for (unsigned j : B) {
                R[j] = reward;
            }
        }

        /* Exploration vector */
        const std::vector<double> Y = generateStochasticVector();

        /* Update row M[island] */
        for (unsigned j = 0; j < nIslands; ++j) {
            migrationMatrix[island][j] = (1.0 - beta) * (alpha * migrationMatrix[island][j] + (1.0 - alpha) * R[j]) + beta * Y[j];
        }

        /* Normalize because of the floating point error */
        normalizeRow(island);
    }

    /* Normalize one row of M */
    void normalizeRow(unsigned island) {
        double sum = 0.0;

        for (double value : migrationMatrix[island]) {
            sum += value;
        }

        if (sum <= 0.0) {
            const double p = 1.0 / static_cast<double>(nIslands);
            for (double& value : migrationMatrix[island]) {
                value = p;
            }
            return;
        }

        for (double& value : migrationMatrix[island]) {
            value /= sum;
        }
    }

    /* Select */
    unsigned selectDestination(unsigned island) const {
        const double random = rndF();
        double cumulative = 0.0;

        for (unsigned j = 0; j < nIslands; ++j) {
            cumulative += migrationMatrix[island][j];

            if (random <= cumulative) {
                return j;
            }
        }

        /* Floating point fallback */
        return nIslands - 1;
    }

    /* Select outgoing particles */
    void selectMigrants(unsigned island, const Swarm& swarm) {
        /* Clear old outgoing packets */
        for (unsigned destination = 0; destination < nIslands; ++destination) {
            outgoing[island][destination].clear();
        }

        /* Every particle participates in the migration decision */
        for (unsigned p = 0; p < swarm.population.size(); ++p) {
            const unsigned destination = selectDestination(island);

            Migrant migrant;
            migrant.particle = swarm.population[p];
            migrant.personalBest = swarm.bestIndividualExp[p];
            migrant.origin = island;
            migrant.fitnessBefore = swarm.population[p].fitness;

            outgoing[island][destination].push_back(std::move(migrant));
        }
    }

    /* Return all migrants destined for an island */
    std::vector<Migrant> collectIncoming(unsigned destination) {
        std::vector<Migrant> incoming;

        for (unsigned source = 0; source < nIslands; ++source) {
            auto& sourcePackets = outgoing[source][destination];
            incoming.insert(incoming.end(), std::make_move_iterator(sourcePackets.begin()), std::make_move_iterator(sourcePackets.end()));
            sourcePackets.clear();
        }

        return incoming;
    }

    /* Print migration matrix */
    void printMatrix(unsigned generation) {
        std::lock_guard<std::mutex> lock(consoleMutex);

        std::printf("\n\n[DIMPSO] Migration matrix after generation %u\n", generation);
        std::printf("             ");

        for (unsigned j = 0; j < nIslands; ++j) {
            std::printf("I%-8u", j);
        }

        std::printf("\n");

        for (unsigned i = 0; i < nIslands; ++i) {
            std::printf("I%-11u", i);
            for (unsigned j = 0; j < nIslands; ++j) {
                std::printf("%.5f ", migrationMatrix[i][j]);
            }
            std::printf("\n");
        }
    }
};

/* Apply incoming migration to an island */
void receiveMigrants(unsigned island, Swarm& swarm, DIMPSOController& controller, ParticleTracking& tracking) {
    std::vector<Migrant> incoming = controller.collectIncoming(island);

    if (incoming.empty()) {
        return;
    }

    /* Construct new population */
    std::vector<Particle> newPopulation;
    std::vector<Particle> newPersonalBest;
    std::vector<unsigned> newOrigin;
    std::vector<double> newFitnessBefore;

    newPopulation.reserve(incoming.size());
    newPersonalBest.reserve(incoming.size());
    newOrigin.reserve(incoming.size());
    newFitnessBefore.reserve(incoming.size());

    for (Migrant& migrant : incoming) {
        newPopulation.push_back(std::move(migrant.particle));
        newPersonalBest.push_back(std::move(migrant.personalBest));
        newOrigin.push_back(migrant.origin);
        newFitnessBefore.push_back(migrant.fitnessBefore);
    }

    /* Replace the island's population */
    swarm.population = std::move(newPopulation);

    swarm.bestIndividualExp = std::move(newPersonalBest);

    /* The social best values are no longer valid beacuse the composition of the swarm changed */
    /* Recreate with the new population size */
    swarm.bestSocialExp.assign(swarm.population.size(), Particle(swarm.nAllele));

    /* Update population size */
    swarm.tPop = static_cast<unsigned>(swarm.population.size());

    /* Update tracking */
    tracking.origin = std::move(newOrigin);

    tracking.fitnessBefore = std::move(newFitnessBefore);
}

/* Run one island */
void runInstance(unsigned islandIndex, Swarm& swarm, unsigned migrationInterval, double alpha, double beta, Barrier& barrier, DIMPSOController& dimpso) {
    /* Unique output names */
    swarm.nfGen = "csvs/" + swarm.nfGen + "_island" + std::to_string(islandIndex);
    swarm.nfRun = "csvs/" + swarm.nfRun;

    const std::string suffix = "_island" + std::to_string(islandIndex);

    if (swarm.nfRun.size() >= 4) {
        swarm.nfRun.insert(swarm.nfRun.size() - 4, suffix);
    } else {
        swarm.nfRun += suffix;
    }

    /* Allocate PSO Memory */
    swarm.initVariables();
    swarm.reserveMemory();

    ParticleTracking tracking;

    /* Independent runs */
    for (unsigned run = 0; run < swarm.nRun; ++run) {
        /* Synchronize all islands */
        barrier.wait();

        /* Controller reset only once */
        if (islandIndex == 0) {
            dimpso.initialize(dimpso.nIslands, alpha, beta);
        }

        barrier.wait();

        /* Restore original population size for every independent run */
        swarm.tPop = static_cast<unsigned>(swarm.population.size());

        /* Statistics */
        initStatistics(swarm.Run);
        
        /* Initial random population */
        swarm.initPopulation();

        /* Reinitialize the personal/social experiences */
        swarm.bestIndividualExp.assign(swarm.tPop, Particle(swarm.nAllele));

        swarm.bestSocialExp.assign(swarm.tPop, Particle(swarm.nAllele));

        /* Initialize tracking */
        tracking.initialize(islandIndex, swarm);

        /* Per-run output. */
        std::string fileGen;

        swarm.runFileName(run, fileGen);

        const unsigned seed = initRandom(0);

        runHeader(fileGen, seed, swarm);

        {
            std::lock_guard<std::mutex> lock(consoleMutex);
            std::printf("\n[DIMPSO ISLAND %u] Run %u started (seed=%u)\n", islandIndex, run, seed);
        }

        /* Main DIMPSO loop */
        for (unsigned gen = 0; gen < swarm.nGen; ++gen) {
            /* STEP 1: Evaluate current population */
            initStatistics(swarm.Gen);
            swarm.evaluatePopulation(gen);

            /* Record statistics */
            swarm.runInfo(fileGen, gen);

            /* STEP 2: Apply the existing PSO algorithm */
            swarm.PSOAlgorithm(gen);

            /* mutation */
            swarm.mutation();

            /* STEP 3: Check migration interval */
            const bool migrationTime = ((gen + 1) % migrationInterval == 0);

            if (!migrationTime) {
                continue;
            }

            /* DIMPSO Migration Epoch */
            /* Analyze current particles */
            initStatistics(swarm.Gen);
            swarm.evaluatePopulation(gen);

            if (islandIndex == 0) {
                dimpso.clearFeedback();
            }

            barrier.wait();

            /* Analyze */
            dimpso.analyze(islandIndex, swarm, tracking);

            barrier.wait();

            /* Learn */
            dimpso.learn();
            
            barrier.wait();

            /* Update */
            dimpso.updateMigrationMatrix(islandIndex);
            
            barrier.wait();

            /* Select */
            dimpso.selectMigrants(islandIndex, swarm);

            barrier.wait();

            /* Migrate */
            receiveMigrants(islandIndex, swarm, dimpso, tracking);

            barrier.wait();

            /* Display M after migration */
            if (islandIndex == 0) {
                dimpso.printMatrix(gen + 1);
            }

            barrier.wait();

            /* Report new population size */
            {
                std::lock_guard<std::mutex> lock(consoleMutex);

                std::printf("\n[DIMPSO] Island %u | Generation %u | Population = %u", islandIndex, gen + 1, swarm.tPop);
            }
        }

        /* Finish Run */
        runFooter(fileGen, swarm);

        {
            std::lock_guard<std::mutex> lock(consoleMutex);

            std::printf("\n[DIMPSO ISLAND %u] Run %u finished | Best fitness = %.3f\n", islandIndex, run, swarm.Run.best.fitness);
        }

        /* Make sure all islands finish this run before beginning another one */
        barrier.wait();
    }

    swarm.freeMemory();
}

} // namespace PSwarm

/* MAIN ENTRY POINT */
int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::printf("Usage: %s <input_file> [num_islands] [migration_interval] [alpha] [beta]\n", argv[0]);
        return 0;
    }

    std::filesystem::create_directories("csvs");

    const std::string inputFile = argv[1];

    /* n Number of islands */
    unsigned numIslands = (argc >= 3) ? static_cast<unsigned>(std::stoul(argv[2])) : 5;

    /* μ Migration Interval */
    unsigned migrationInterval = (argc >= 4) ? static_cast<unsigned>(std::stoul(argv[3])) : 100;

    /* α Alpha parameter */
    double alpha = (argc >= 5) ? std::stod(argv[4]) : 0.8;

    /* β Beta parameter */
    double beta = (argc >= 6) ? std::stod(argv[5]) : 0.01;

    /* Validate parameters. */
    if (numIslands < 2) {
        std::printf("Error: num_islands must be >= 2.\n");
        return 1;
    }

    if (migrationInterval == 0) {
        std::printf("Error: migration_interval must be > 0.\n");
        return 1;
    }

    if (alpha <= 0.0 || alpha >= 1.0) {
        std::printf("Error: alpha must be in the range [0, 1].\n");
        return 1;
    }

    if (beta <= 0.0 || beta >= 1.0) {
        std::printf("Error: beta must be in the range [0, 1].\n");
        return 1;
    }

    std::printf(
        "\n"
        "====================================================\n"
        "Dynamic Island Model Particle Swarm Optimization\n"
        "====================================================\n"
        "Input file            : %s\n"
        "Number of islands     : %u\n"
        "Migration interval μ  : %u\n"
        "Alpha α               : %.4f\n"
        "Beta β                : %.4f\n"
        "Migration Topology    : Fully Connected\n"
        "Migration Policy      : Dynamic matrix M\n"
        "Learning Matrix       : Dm\n"
        "Migration Rate        : Deterined by \n"
        "====================================================\n\n",
        inputFile.c_str(),
        numIslands,
        migrationInterval,
        alpha,
        beta
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

    /* DIMPSO Controller */
    PSwarm::DIMPSOController dimpso(numIslands, alpha, beta);

    /* Synchronization Barrier */
    PSwarm::Barrier migrationBarrier(numIslands);

    /* Launch n Islands */
    std::vector<std::thread> workers;
    workers.reserve(numIslands);

    for (unsigned i = 0; i < numIslands; ++i) {
        workers.emplace_back(PSwarm::runInstance, i, std::ref(*islands[i]), migrationInterval, alpha, beta, std::ref(migrationBarrier), std::ref(dimpso));
    }

    /* Wait for all sub-swarms to finish */
    for (auto& worker : workers) {
        worker.join();
    }

    std::printf("\nAll %u DIMPSO islands finished.\n", numIslands);

    return 0;
}