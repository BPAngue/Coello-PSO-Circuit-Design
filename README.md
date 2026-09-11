# PSO Matrix Circuit

A Particle Swarm Optimization (PSO) implementation for circuit matrix design, based on Carlos Coello's original work. This repository includes the original Spanish-language source, an English translation, and a C++ rendition.

## Attribution

This project is based on the Particle Swarm Optimization (PSO) circuit design algorithm originally developed by Carlos A. Coello Coello, Erika Luna, and Arturo Hernandez-Aguirre. The `orig_code/` folder preserves the original Spanish-language implementation; `translated/` and `cplusplus/` are derivative adaptations for readability and portability. Please cite the original work if you build on or publish results using this code:

```bibtex
@inbook{coello2003pso,
  author    = {Coello, Carlos and Luna, Erika and Hernandez-Aguirre, Arturo},
  title     = {Use of Particle Swarm Optimization to Design Combinational Logic Circuits},
  journal   = {Lecture Notes in Computer Science},
  volume    = {2606},
  pages     = {398-409},
  year      = {2003},
  month     = {06},
  isbn      = {978-3-540-00730-2},
  doi       = {10.1007/3-540-36553-2_36}
}
```

## Repository Structure

| Folder | Description |
|---|---|
| `orig_code/` | Original Coello code (Spanish variable/function names, `circuitos.c`, `estadisticas.c`) |
| `translated/` | English translation of the original Coello code (`circuits.c`, `statistics.c`) |
| `cplusplus/` | C++ renditions of the codebase (more details of the folder structure is include inside the folder) |
| `data_files/` | Input data files (`.dta`) used with the original Coello code |

### Translated Code (English)

```bash
gcc -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuits.c statistics.c random.c -lm
```

### Original Coello Code (Spanish)

```bash
gcc -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuitos.c estadisticas.c random.c -lm
```

### The Earlier Rendition of the C++ Variant (no Makefile) (NOT CURRENTLY APPLICABLE)

Standard build:

```bash
g++ -O2 -o psomatrixcircuit.exe psomatrixcircuit.cpp matrixpso.cpp circuits.cpp statistics.cpp random.cpp
```

Optimized build (release mode):

```bash
g++ -std=c++17 -O3 -DNDEBUG -march=native -ffast-math -o psomatrixcircuit.exe psomatrixcircuit.cpp matrixpso.cpp circuits.cpp statistics.cpp random.cpp
```

### Troubleshooting: Linker Errors

If you encounter common symbol / multiple-definition linker errors in the C builds, add the `-fcommon` flag:

**Translated code:**

```bash
gcc -fcommon -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuits.c statistics.c random.c -lm
```

**Original Coello code:**

```bash
gcc -fcommon -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuitos.c estadisticas.c random.c -lm
```

## Usage

Run the compiled program with an input data file:

```bash
./psomatrixcircuit inputfile.dta

For `cplusplus/Original/`, use the `.dta` files under `cplusplus/Original/data_files/` (packed PLU truth table format).