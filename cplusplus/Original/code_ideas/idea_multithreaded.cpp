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
        }
        else {
            cv.wait(lock, [&]() {
                return generation != currentGeneration;
            });
        }
    }
};

/*
 * ============================================================
 * Shared state between all sub-swarms
 * ============================================================
 */
struct ParallelState {
    
    unsigned numSwarms;

    /*
     * Each element contains the Gen.best of one sub-swarm.
     * 
     * swarmBests[0] = Swarm 0 Gen.best
     * swarmBests[1] = Swarm 1 Gen.best
     * ...
     */
    std::vector<Particle> swarmBests;

    /*
     * Best Gen.best among all sub-swarms
     */
    Particle globalBest;

    /*
     * Synchronization point
     */
    Barrier barrier;

    /*
     * Protects access to swarmBests/globalBest while they
     * are being written/read.
     */
    std::mutex bestMutex;

    ParallelState(unsigned count) : numSwarms(count), swarmBests(count), barrier(count) {}
};

/*
 * ============================================================
 * Find the global best solution among all swarmBests
 * ============================================================
 */
void calculateSwarmBest(unsigned instanceIndex, Swarm& swarm, ParallelState& state)
{
    /*
     * Store this sub-swarm's best solution.
     */
    {
        std::lock_guard<std::mutex> lock(state.bestMutex);

        state.swarmBests[instanceIndex] = swarm.Gen.best;
    }

    /*
     * Wait until EVERY sub-swarm has written its Gen.best.
     */
    state.barrier.wait();

    /*
     * Only one thread needs to calculate the best.
     * 
     * We use instance 0 as the coordinator.
     */
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

    /*
     * Wait until instance 0 has finished calculating globalBest
     */
    state.barrier.wait();
}

/*
 * ============================================================
 * Apply swarm-best to every sub-swarm
 * ============================================================
 */
void applyGlobalBest(Swarm& swarm, ParallelState& state)
{
    const double phi3 = 1.0;  // hardcoded global-best coefficient

    for (unsigned i = 0; i < swarm.tPop; ++i) {
        for (unsigned d = 0; d < swarm.nAllele; ++d) {

            const double phi1p = rndF() * swarm.phi1;
            const double phi2p = rndF() * swarm.phi2;
            const double phi3p = rndF() * phi3;

            /*
             * Single-evolutionary phase:
             *
             * pbest = bestIndividualExp[i]
             * sbest = bestSocialExp[i]
             * gbest = state.globalBest
             */
            swarm.population[i].vi[d] +=
                phi1p *
                (
                    static_cast<double>(
                        swarm.bestIndividualExp[i].chromX[d]
                    )
                    -
                    static_cast<double>(
                        swarm.population[i].chromX[d]
                    )
                );

            swarm.population[i].vi[d] +=
                phi2p *
                (
                    static_cast<double>(
                        swarm.bestSocialExp[i].chromX[d]
                    )
                    -
                    static_cast<double>(
                        swarm.population[i].chromX[d]
                    )
                );

            swarm.population[i].vi[d] +=
                phi3p *
                (
                    static_cast<double>(
                        state.globalBest.chromX[d]
                    )
                    -
                    static_cast<double>(
                        swarm.population[i].chromX[d]
                    )
                );

            /* Clamp velocity */
            swarm.population[i].vi[d] =
                std::clamp(
                    swarm.population[i].vi[d],
                    -swarm.vMax,
                    swarm.vMax
                );

            /* Update position via sigmoid-mapped velocity */
            const double vNorm =
                sigmoid(swarm.population[i].vi[d]);

            switch (swarm.representation) {

                case BINARY:
                    swarm.population[i].chromX[d] =
                        flip(vNorm);
                    break;

                case INTEGER_A:
                    swarm.population[i].chromX[d] =
                        flip(vNorm)
                            ? swarm.bestSocialExp[i].chromX[d]
                            : swarm.population[i].chromX[d];
                    break;

                case INTEGER_B:
                    swarm.population[i].chromX[d] =
                        flip(vNorm)
                            ? swarm.bestSocialExp[i].chromX[d]
                            : flip(1.0 - vNorm)
                                ? swarm.bestIndividualExp[i].chromX[d]
                                : swarm.population[i].chromX[d];
                    break;
            }
        }
    }
}

/*
 * ============================================================
 * Run one sub-swarm
 * ============================================================
 */
void runInstance(unsigned instanceIndex, Swarm& swarm, ParallelState& state, unsigned M, unsigned N) 
{
    /*
     * Unique output filenames.
     */
    swarm.nfGen = "csvs/" + swarm.nGen;
    swarm.nfGen += "_swarm" + std::to_string(instanceIndex) + "_";

    swarm.nfRun = "csvs/" + swarm.nfRun;

    const std::string suffix = "_swarm" + std::to_string(instanceIndex);

    if (swarm.nfRun.size() >= 4) {
        swarm.nfRun.insert(
            swarm.nfRun.size() - 4,
            suffix
        );
    }

    /*
     * Allocate PSO memory
     */
    swarm.initVariables();
    swarm.reserveMemory();

    /*
     * One thread owns one Swarm object.
     */
    globalHeader(swarm.nfRun, swarm);

    /*
     * ========================================================
     * Independent runs
     * ========================================================
     */
    for (unsigned run = 0; run < swarm.nRun; ++run) {
        {
            std::lock_guard<std::mutex> lock(consoleMutex);

            std::printf("\n[Swarm %u] Run %02u Started\n", instanceIndex, run);
        }

        /*
         * ----------------------------------------------------
         * Initialize this run
         * ----------------------------------------------------
         */
        initStatistics(swarm.Run);
        swarm.initPopulation();

        /*
         * ====================================================
         * OUTER LOOP
         *
         * This corresponds to:
         *
         *          M iterations?
         *
         * in your paper's flowchart.
         * ====================================================
         */
        for (unsigned m = 0; m < M; ++m) {
            {
                std::lock_guard<std::mutex> lock(consoleMutex);

                if (instanceIndex == 0) {
                    std::printf("\n========== Outer iteration %u / %u ==========\n", m + 1, M);
                }
            }

            /*
             * =================================================
             * LOCAL PSO EVOLUTION
             *
             * Each sub-swarm evolves independently for N
             * iterations.
             * =================================================
             */
            for (unsigned n = 0; n < N; ++n) {
                /*
                 * IMPORTANT:
                 *
                 * Do NOT use n as the generation argument.
                 *
                 * evaluatePopulation(0) means:
                 *
                 *     "initialize/update personal best"
                 *
                 * and would reset the meaning of personal-best
                 * at every outer iteration.
                 *
                 * We therefore use a monotonically increasing
                 * generation number.
                 */
                const unsigned generation = m * N + n;

                initStatistics(swarm.Gen);
                swarm.evaluatePopulation(generation);
                
                /*
                 * Update normal run statistics.
                 *
                 * We can use an empty filename here if we do
                 * not want generation CSV output from this
                 * experimental version.
                 */
                swarm.runInfo("", generation);

                /*
                 * Normal LOCAL PSO update:
                 *
                 * Personal best
                 *       +
                 * Local neighborhood best
                 */
                swarm.PSOAlgorithm(generation);

                /*
                 * Existing mutation remains enabled.
                 */
                swarm.mutation();
            }

            /*
             * =================================================
             * SYNCHRONIZATION POINT
             * =================================================
             *
             * Every sub-swarm has now completed N local
             * iterations.
             */
            calculateSwarmBest(instanceIndex, swarm, state);

            /*
             * =================================================
             * SWARM-BEST PHASE
             * =================================================
             */
            applyGlobalBest(swarm, state);

            /*
             * =================================================
             * SINGLE EVOLUTIONARY PHASE
             *
             * The existing PSOAlgorithm() now uses:
             *
             *     personal best
             *          +
             *     swarm best
             *
             * instead of:
             *
             *     personal best
             *          +
             *     local neighborhood best
             * =================================================
             */
            swarm.PSOAlgorithm(
                m * N + N
            );


            /*
             * Mutation after the swarm-best PSO update.
             */
            // swarm.mutation();


            /*
             * Wait so that every swarm finishes the
             * swarm-best phase before any swarm begins
             * the next outer iteration.
             */
            state.barrier.wait();
        }

        /*
         * ====================================================
         * End of run
         * ====================================================
         */

        {
            std::lock_guard<std::mutex> lock(consoleMutex);

            std::printf("\n[Swarm %u] Run %02u Finished\n", instanceIndex, run);
        }
    }

    swarm.freeMemory();

    {
        std::lock_guard<std::mutex> lock(consoleMutex);

        std::printf(
            "[Swarm %u] Complete.\n",
            instanceIndex
        );
    }
} 

} // namespace PSwarm

/*
 * ============================================================
 * MAIN
 * ============================================================
 */
int main(int argc, char* argv[])
{
    if (argc < 2) {

        std::printf(
            "Usage: %s <input_file> [num_swarms] [M] [N]\n",
            argv[0]
        );

        return 0;
    }


    std::filesystem::create_directories("csvs");


    const std::string inputFile = argv[1];


    /*
     * K = number of parallel sub-swarms.
     */
    unsigned numSwarms = (argc >= 3) ? static_cast<unsigned>(std::stoul(argv[2])) : std::thread::hardware_concurrency();

    if (numSwarms == 0)
        numSwarms = 1;

    /*
     * M = number of outer parallel-swarm iterations.
     *
     * Default = 10.
     */
    unsigned M =
        (argc >= 4)
        ? static_cast<unsigned>(std::stoul(argv[3]))
        : 20;


    /*
     * N = number of independent local PSO generations
     * performed by each sub-swarm before synchronization.
     *
     * Default = 5.
     */
    unsigned N =
        (argc >= 5)
        ? static_cast<unsigned>(std::stoul(argv[4]))
        : 10;


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


    /*
     * --------------------------------------------------------
     * Create K independent Swarm objects.
     *
     * Each object represents one sub-swarm.
     * --------------------------------------------------------
     */
    std::vector<std::unique_ptr<PSwarm::Swarm>> swarms;

    swarms.reserve(numSwarms);


    for (unsigned i = 0; i < numSwarms; ++i)
    {
        auto swarm = std::make_unique<PSwarm::Swarm>();

        if (!swarm->loadParameters(inputFile)) {

            std::printf(
                "Failed to load parameters for swarm %u.\n",
                i
            );

            return 1;
        }

        swarms.push_back(std::move(swarm));
    }


    /*
     * Shared synchronization / communication state.
     */
    PSwarm::ParallelState state(numSwarms);


    /*
     * --------------------------------------------------------
     * Launch K sub-swarms.
     * --------------------------------------------------------
     */
    std::vector<std::thread> workers;

    workers.reserve(numSwarms);


    for (unsigned i = 0; i < numSwarms; ++i)
    {
        workers.emplace_back(
            PSwarm::runInstance,
            i,
            std::ref(*swarms[i]),
            std::ref(state),
            M,
            N
        );
    }


    /*
     * --------------------------------------------------------
     * Wait for all sub-swarms.
     * --------------------------------------------------------
     */
    for (auto& worker : workers)
        worker.join();


    std::printf(
        "\nAll %u sub-swarms finished.\n",
        numSwarms
    );


    return 0;
}