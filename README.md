# Compilation Process
Compile the program with all of the needed header files:

## For Translated Code

```gcc -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuits.c statistics.c random.c -lm```

## For Original Coello Code

```gcc -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuitos.c estadisticas.c random.c -lm```

If the previous line does not work, use the following line instead:

## For Translated Code

```gcc -fcommon -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuits.c statistics.c random.c -lm```

## For Original Coello Code

```gcc -fcommon -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuitos.c estadisticas.c random.c -lm```

Run the program with the input data

```./psomatrixcircuit inputfile.dta```