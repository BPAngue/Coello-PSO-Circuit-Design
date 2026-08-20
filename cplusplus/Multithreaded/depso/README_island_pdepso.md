# Island-Model PDEPSO with Particle Speed Limit

Implements the pipeline from:

> Ikegami, H. and Mori, H. (2018). "Development of DEPSO Island Model with
> Particle Speed Limit for Distribution Network Reconfigurations."
> IFAC PapersOnLine 51-28, pp. 552-557.

applied to your existing gate-matrix circuit-design PSO codebase.

## New files (drop-in additions, nothing else changes)

- `epso_island.h` / `epso_island.cpp` — the algorithm itself:
  - **EPSO** (eqns 8–11/16–19): each particle's acceleration weights `w1`,
    `w2` are mutated with Gaussian noise each generation; `replicationRate`
    mutated children are generated per particle and the best of
    `{parent, children}` is kept ("(µ+λ)" selection, so fitness never
    regresses).
  - **Particle speed limit** (eqn 22): velocity is clamped to
    `[-vMax, +vMax]` before the sigmoid decode, reusing `swarm.vMax`.
  - **Island model** (Fig. 3–4, section 4.2): the population is split into
    `nIslands` independent sub-populations that evolve separately and
    exchange their best particles in a ring topology every
    `migrationInterval` generations, replacing the worst particles of the
    downstream island (`migrationRate` controls how many).
- `run_island_pdepso.cpp` — a new `main()` that loads a circuit data file
  (**same file format** your existing tool uses) and runs the island
  pipeline instead of the original single-population PSO.
- `psomatrixcircuit.cpp` — **one small change**: the original `main()` is
  now wrapped in `#ifndef PSO_NO_MAIN` so this file's non-`main` code
  (`Swarm::loadParameters`, `initVariables`, `evaluateParticle`, etc.) can
  be linked into both the original executable and the new island demo
  without a duplicate-`main` conflict. Nothing else in the file changed.

Everything else (`random.*`, `circuits.*`, `matrixpso.*`, `statistics.*`,
`particle.h`, `psomatrixcircuit.h`) is untouched.

## Design decisions worth knowing about

- **Representation requirement**: the paper encodes everything as a flat
  bit-string decoded via the sigmoid function, matching your existing
  `BINARY` representation. `IslandPDEPSO::execute()` checks
  `swarm.representation == BINARY` and prints an error otherwise — the
  `INTEGER_A`/`INTEGER_B` representations aren't part of this paper's
  formulation.
- **gbest mutation on a discrete chromosome**: eqn (19),
  `g*best = gbest + τ'·N(0,1)`, perturbs a *continuous* gbest vector in the
  original EPSO. Since your gbest is a discrete chromosome, that term is
  instead realized as an extra additive noise contribution in the velocity
  update (see comments in `epso_island.h`/`.cpp`) — it plays the same
  "randomize the pull toward gbest" role without corrupting the discrete
  chromosome bits directly.
- **Output format is unchanged**: `IslandPDEPSO::execute()` reuses your
  existing `globalHeader`/`runHeader`/`generationStats` (via
  `Swarm::runInfo`)/`runFooter`/`runStatistics` free functions, so the CSV
  files it produces are identical in format to the ones from the original
  `psomatrixcircuit` executable — just with the different search dynamics.

## Building

```bash
# Library objects (unchanged files)
g++ -std=c++17 -O2 -c random.cpp circuits.cpp matrixpso.cpp statistics.cpp

# psomatrixcircuit.cpp's non-main code, built for reuse in the island demo
g++ -std=c++17 -O2 -DPSO_NO_MAIN -c psomatrixcircuit.cpp -o psomatrixcircuit_lib.o

# New island-PDEPSO executable
g++ -std=c++17 -O2 \
    random.o circuits.o matrixpso.o statistics.o psomatrixcircuit_lib.o \
    epso_island.cpp run_island_pdepso.cpp -o pdepso_island

# (Original executable still builds exactly as before, main included)
g++ -std=c++17 -O2 random.o circuits.o matrixpso.o statistics.o \
    psomatrixcircuit.cpp -o psomatrixcircuit
```

## Running

```bash
./pdepso_island <input_file> [nIslands migInterval migRate replication tau tauPrime]
```

All six tuning arguments are optional; omitted ones fall back to the
paper's Table 1 defaults for Methods C/F:

| Parameter          | Default | Meaning                                   |
|--------------------|---------|--------------------------------------------|
| `nIslands`         | 4       | number of sub-populations                  |
| `migrationInterval`| 400     | generations between migrations             |
| `migrationRate`    | 0.1     | fraction of each island migrated           |
| `replicationRate`  | 2       | EPSO children generated per particle       |
| `tau`              | 0.01    | learning rate for w1/w2 mutation           |
| `tauPrime`         | 0.02    | learning rate for gbest-noise mutation     |

The input file format is **unchanged** — the same `.txt` data file you use
with the original `psomatrixcircuit` executable works here too, as long as
`representation` is set to `2` (Binary).

This was verified end-to-end on a small synthetic 2-input XOR circuit: it
converges to a 0-violation solution and produces CSV output in the same
format as the original tool.
