# Compilation Process
Compile the program with all of the needed header files:

```gcc -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuits.c statistics.c random.c -lm```

If the previous line does not work, use the following line instead:

```gcc -fcommon -o psomatrixcircuit psomatrixcircuit.c matrixpso.c circuits.c statistics.c random.c -lm```

Run the program with the input data

```./psomatrixcircuit inputfile.dta```