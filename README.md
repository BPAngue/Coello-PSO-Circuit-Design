# Compilation Process
Compile the program with all of the needed header files:

## For Translated Code

```gcc -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuits.c statistics.c random.c -lm```

## For Original Coello Code

```gcc -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuitos.c estadisticas.c random.c -lm```

## For C++ Rendition of the Code

```g++ -O2 -o psomatrixcircuit.exe psomatrixcircuit.cpp matrixpso.cpp circuits.cpp statistics.cpp random.cpp```

```g++ -std=c++17 -O3 -DNDEBUG -march=native -ffast-math -o psomatrixcircuit.exe psomatrixcircuit.cpp matrixpso.cpp circuits.cpp statistics.cpp random.cpp```

If you encounter common linker errors, add -fcommon flag.

## For Translated Code

```gcc -fcommon -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuits.c statistics.c random.c -lm```

## For Original Coello Code

```gcc -fcommon -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuitos.c estadisticas.c random.c -lm```

Run the program with the input data

```./psomatrixcircuit inputfile.dta```