# README #

InterLig is a software for the sequence-order alignment and Virtual Screening of small molecules.


## Installation ##

To compile InterLig just clone this repo and run the command:

```
make
```

A binary named InterLig will be generated in the cloned folder.

## Usage ##
To run the binary in the folder where it is installed:

    ./InterLig -mol2 target.mol2 database.mol2 <options>

    ### Main options: ###
    
    -mol2  : Compare a mol2 file with a mol2 database (multiple molecules in the database)
             (./InterLig -mol2 target.mol2 database.mol2)

    ### Other options: ###
    
    -d0           : d0 parameter in the Levitt-Gerstein score (default: 0.5)
    -dW           : weight of sequence similarity (e.g. BLOSUM62) in the optimization procedure (default: 0.5)
    -eps          : weight of structural similarity in the optimization procedure (default: 0.5)
    -nullP        : percentage of ignored atoms in the smaller molecule (default: 0)
    -v            : verbose output
    -anneal <int> : number of annealing rounds (default: 1)
    -seed <int>   : seed for the stochastic process (random number generation)
    -ch <int>     : Markov chain length multiplier (default: 10)\n
    -matrix <path>: optional path to the sequence similarity matrix (e.g. BLOSUM62)
    -super <path> : calculate optimal superposition between two aligned PDB files (only -PDB mode)
                    (output to <patg>)
    -h            : print this help

## Usage examples ##

If we want to run on mol2 files with d0 = 0.4 and eps = 0.75:

```
./InterLig -mol2 examples/target.mol2 examples/database.mol2 -d0 0.4 -eps 0.75
```

The database file here is *always* the second file following the mol2 flag.
