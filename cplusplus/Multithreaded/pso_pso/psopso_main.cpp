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
/* File: psopso_main.cpp                                */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Entry point for the Parallel Swarms      */
/* Oriented PSO (PSO-PSO) pipeline. Reuses the same input */
/* data file format as psomatrixcircuit / psomatrixcircuit_mt */
/* (truth table + PSO parameters); the PSO-PSO-specific    */
/* topology parameters (K, N, M, cycles, phi3) are given   */
/* on the command line rather than added to that file, so  */
/* existing input files keep working unmodified with every */
/* program in this project.                                */
/*                                                          */
/* Usage:                                                   */
/*   psopso <input_file> [K] [N] [M] [cycles] [phi3]        */
/*                                                            */
/*   K       number of sub-swarms                (default 4)  */
/*   N       multi-evolutionary phase iterations (default 20)  */
/*   M       single-evolutionary phase iterations(default 10)   */
/*   cycles  number of (multi,single) alternations (default 10)  */
/*   phi3    acceleration toward gbest, c3        (default 2.0)   */
/********************************************************/

#include "psopso.h"

#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::printf(
            "Usage: %s <input_file> [K] [N] [M] [cycles] [phi3]\n\n"
            "  Implements the Parallel Swarms Oriented PSO (PSO-PSO)\n"
            "  algorithm of Gonsalves & Egashira (2013) on top of the\n"
            "  circuit-design Swarm.\n\n"
            "  K       number of sub-swarms                (default 4)\n"
            "  N       multi-evolutionary phase iterations (default 20)\n"
            "  M       single-evolutionary phase iterations (default 10)\n"
            "  cycles  number of (multi,single) alternations (default 10)\n"
            "  phi3    acceleration toward gbest, c3         (default 2.0)\n",
            argv[0]);
        return 0;
    }

    PSwarm::PSOPSOConfig cfg;
    cfg.inputFile        = argv[1];
    cfg.numSwarms        = argc > 2 ? static_cast<unsigned>(std::stoul(argv[2])) : 4;
    cfg.multiPhaseIters  = argc > 3 ? static_cast<unsigned>(std::stoul(argv[3])) : 20;
    cfg.singlePhaseIters = argc > 4 ? static_cast<unsigned>(std::stoul(argv[4])) : 10;
    cfg.numCycles        = argc > 5 ? static_cast<unsigned>(std::stoul(argv[5])) : 10;
    cfg.phi3             = argc > 6 ? std::stod(argv[6]) : 2.0;

    if (cfg.numSwarms == 0) cfg.numSwarms = 1;

    PSwarm::runPSOPSO(cfg);
    return 0;
}
