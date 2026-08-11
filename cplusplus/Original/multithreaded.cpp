#include "psomatrixcircuit.h"
#include "random.h"

#include <string>
#include <thread>
#include <vector>
#include <mutex>
#include <atomic>
#include <filesystem>

namespace PSwarm {

std::mutex consoleMutex;

void runInstance(unsigned instanceIndex, const std::string& inputFile) {
    
    PSwarm::Swarm swarm;

    if (!swarm.loadParameters(inputFile)) {
        std::lock_guard<std::mutex> lock(consoleMutex);
        std::printf("[instance %u] failed to load parameters, aborting.\n", instanceIndex);
        return;
    }

    /* Make unique output files for each instance */
    swarm.nfGen = "csvs/" + swarm.nfGen;
    swarm.nfGen += "_inst" + std::to_string(instanceIndex) + "_";

    swarm.nfRun = "csvs/" + swarm.nfRun;
    swarm.nfRun.insert(swarm.nfRun.size() - 4, "_inst" + std::to_string(instanceIndex));

    swarm.initVariables();
    swarm.reserveMemory();

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

} // namespace PSwarm


int main(int argc, char *argv[]) {

    if (argc < 2) {
        std::printf("Usage: %s <input_file> [num_instances]\n", argv[0]);
        return 0;
    }

    std::filesystem::create_directories("csvs");

    const std::string inputFile = argv[1];

    // check if arguments contain number of threads, if not use the device's max number of threads
    unsigned numInstances = (argc >= 3) ? static_cast<unsigned>(std::stoul(argv[2])) : std::thread::hardware_concurrency();

    // make sure that there is at least 1 thread working.
    if (numInstances == 0) {
        numInstances = 1;
    }

    std::printf("Launching %u independent solver instance(s) on '%s'\n", numInstances, inputFile.c_str());

    std::vector<std::thread> workers;
    workers.reserve(numInstances);

    for (unsigned i = 0; i < numInstances; ++i) {
        workers.emplace_back(PSwarm::runInstance, i, inputFile);
    }

    for (auto& t : workers) {
        t.join();
    }

    std::printf("All %u instance(s) finished.\n", numInstances);

    return 0;
}