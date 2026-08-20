/********************************************************/
/* File: run_island_pdepso.cpp                           */
/* Description: Entry point that loads a circuit data     */
/* file exactly like the original psomatrixcircuit.cpp    */
/* does, then runs Parallel Discrete EPSO with an island   */
/* model and particle speed limit (Ikegami & Mori, 2018)   */
/* instead of the original single-population PSO.          */
/*                                                        */
/* Usage:                                                 */
/*   pdepso_island <input_file> [nIslands migInterval     */
/*                                migRate replication      */
/*                                tau tauPrime]            */
/*                                                        */
/* The input file format is unchanged -- the same data     */
/* file that works with psomatrixcircuit still works       */
/* here; representation must be Binary (BINARY=2) for      */
/* this pipeline.                                          */
/********************************************************/

#include "psomatrixcircuit.h"
#include "epso_island.h"
#include "random.h"

#include <cstdio>
#include <cstdlib>

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::printf("Usage: %s <input_file> "
                     "[nIslands migInterval migRate replication tau tauPrime]\n",
                     argv[0]);
        return 0;
    }

    PSwarm::Swarm swarm;

    if (!swarm.loadParameters(argv[1])) {
        return 0;
    }

    swarm.initVariables();
    swarm.reserveMemory();

    PSwarm::IslandParams params; /* defaults = Table 1, Methods C/F */

    if (argc > 2) params.nIslands          = static_cast<unsigned>(std::atoi(argv[2]));
    if (argc > 3) params.migrationInterval = static_cast<unsigned>(std::atoi(argv[3]));
    if (argc > 4) params.migrationRate     = std::atof(argv[4]);
    if (argc > 5) params.replicationRate   = static_cast<unsigned>(std::atoi(argv[5]));
    if (argc > 6) params.tau               = std::atof(argv[6]);
    if (argc > 7) params.tauPrime          = std::atof(argv[7]);

    std::printf("\nIsland-PDEPSO parameters:\n");
    std::printf("  Islands:            %u\n", params.nIslands);
    std::printf("  Migration interval: %u\n", params.migrationInterval);
    std::printf("  Migration rate:     %.3f\n", params.migrationRate);
    std::printf("  Replication rate:   %u\n", params.replicationRate);
    std::printf("  tau, tau':          %.4f, %.4f\n", params.tau, params.tauPrime);
    std::printf("  Vmax (speed limit): %.3f\n\n", swarm.vMax);

    PSwarm::IslandPDEPSO pdepso(swarm, params);
    pdepso.execute();

    swarm.freeMemory();
    return 1;
}
