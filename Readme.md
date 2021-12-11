Biorseo (Bi-Objective RNA Structure Efficient Optimizer)
===================================

This tool predicts the secondary structure of a RNA sequence with pieces of 3D information (non-canonical contacts) at some places, 
by identifying zones that can fold like known modules from data like the RNA 3D Motif Atlas or Rna3Dmotifs.

Contact : louis.becquey@univ-evry.fr (deprecated), fariza.tahi@univ-evry.fr (head of team).

1/ How it works
===================================
INPUT:
- An RNA sequence (you can go up to ~230 bases, also depends on your hardware, patience, and module average size)

THEN
- **Pattern-matching step** : Find all possible occurrences of known RNA modules in the query sequence, by finding subsequences of the query that match the modules. Our approach is the regular-expression based approach from RNA-MoIP (*Reinharz et al, 2012*), we deal the same way with short components and wildcards. This step can also be delegated to external tools (like Jar3d or BayesPairing) which will use probabilistic models to score the modules insertion sites.
- **Constraints  definition  step** : Define constraints on the secondary structure imposed by modules if they would be included (example, the loop closing base-pairs for modules that are secondary-structure elements).
- **Solve a bi-objective IP problem** : Find a secondary structure that satisfies as much as possible both an energy-based criterion (e.g. minimize energy, or maximize expected accracy), and a criterion taking into account module inclusions, by solving a bi-objective integer linear programming problem, using the previous constraints defined in the previous step.

OUTPUT:
- A set of secondary structures from the Pareto front,
- The list of known modules inserted inplace in the corresponding structures,
- The scores of the secondary structures on the used objective functions.

2/ The different models
==================================
### MODULE SOURCES

See supported module sources in file [SOURCES.md](SOURCES.md).

### OBJECTIVE FUNCTIONS TO OPTIMIZE

First, the energy-based objective functions :
* **MEA** : Try to maximize the expected accuracy of the secondary structure. This should lead to the MEA secondary structure estimator. Note that in practice, Biorseo maximizes the sum of observed basepairs' probabilities, therefore the scores returned are not interpretable as expected accuracies.
* **MFE** : Try to minimize the free energy of the secondary structure. The free energy model is very simple and considers sequence-dependent contributions for each stack of basepairs. This should lead to the MFE secondary structure estimator. Note that in practice, Biorseo maximizes the opposite of the energy, so the given score can be seen as the opposite of energy (which should be negative).

These approaches are the IPknot (*Sato et al. 2011*) and Biokop (*Legendre et al. 2018*) approaches.

Then, the objective functions for the insertion of modules :
* **Function A** : weights a module by its squared number of nucleotides (like RNA-MoIP).
* **Function B** : weights a module by its number of components (strands) and penalizes it by the log^(_2) of its nucleotide size.
* **Function C** : weights a module by its insertion site score (JAR3D or BayesPairing score).
* **Function D** : weights a module by its number of components (strands) and insertion site score (JAR3D or BayesPairing score), and penalizes it by the log^(_2) of its nucleotide size.


3/ Installation
==================================
You can install Biorseo using Docker, or compile it from source. Check the file [INSTALL.md](INSTALL.md) for installation instructions.


4/ How to run
==================================
The command to run is different depending on your installation method, see [INSTALL.md](INSTALL.md). But both methods allow you to define options to get the desired behavior. 

```
Usage:  You must provide:
        1) a FASTA input file with -s,
        2) a module type with --rna3dmotifs, --carnaval, --json or --pre-placed,
        3) one module-based scoring function with --func A, B, C, or D
        4) one energy-based scoring function with --mfe or --mea,
        5) how to display results: in console (-v), or in a result file (-o).

Options:
  -h [ --help ]                     Print the help message
  --version                         Print the program version
  -s [ --seq ] arg                  Fasta file containing the RNA sequence
  -d [ --descfolder ] arg           A folder containing modules in .desc 
                                    format, as produced by Djelloul & Denise's 
                                    catalog program (deprecated)
  -r [ --rinfolder ] arg            A folder containing CaRNAval's RINs in .txt
                                    format, as produced by script 
                                    transform_caRNAval_pickle.py
  -j [ --jsonfolder ] arg           A folder containing a custom motif library 
                                    in .json format
  -x [ --pre-placed ] arg           A CSV file providing motif insertion sites 
                                    obtained with another tool.
  -f [ --function ] arg (=B)        (A, B, C, or D) Objective function to score
                                    module insertions:
                                      (A) insert big modules
                                      (B) light, high-order modules
                                      (C) well-scored modules
                                      (D) light, high-order, well-scored
                                        modules
                                      C and D require position scores
                                      provided by --pre-placed.
                                    
  -E [ --mfe ]                      Minimize stacking energies
                                      (leads to MFE extimator)
  -A [ --mea ]                      (default) Maximize expected accuracy
                                      (leads to MEA estimator)
  -c [ --first-objective ] arg (=2) (1 or 2) Objective to solve in the 
                                    mono-objective portions of the algorithm.
                                      (1) is the module objective,
                                      given by --function
                                      (2) is energy-based objective,
                                      either MFE or MEA
  -o [ --output ] arg               A file to summarize the computation results
  -t [ --theta ] arg (=0.001)       Pairing probability threshold to consider 
                                    or not the possibility of pairing
  -n [ --disable-pseudoknots ]      Add constraints forbidding the formation of
                                    pseudoknots
  -l [ --limit ] arg (=500)         Intermediate number of solutions in the 
                                    Pareto set above which we give up the 
                                    calculation.

```
## Run in Biokop-mode :
If you run Biorseo with both options `--mfe` and `--mea`, the biobjective optimisation problem will be defined without modules, comparing the two energy-based criterions together. This should be equivalent to run the Biokop tool (*Legendre et al. 2018*) with only one Pareto-set (option `biokop -n 1`).

5/ Example output and interpretation
================================

Let's consider an example input FASTA file (data/fasta/example.fa):
```
>PDB_00376
GGGGUAUCGCCAAGCGGUAAGGCACCGGAUUCUGAUUCCGGAGGUCGAGGUUCGAAUCCUCGUACCCCAGCCA
```

By running `docker run -v `pwd`/data:/workdir/data persalteas/biorseo -s data/fasta/example.fa --descfolder data/modules/DESC --func B -v -o data/MyOutput.biorseo`, i get many information about what is happening:

### BASEPAIR PROBABILITIES
```
Summary of basepair probabilities:

        === -log10(p(i,j)) for each pair (i,j) of nucleotides: ===

        GGGGUAUCGCCAAGCGGUAAGGCACCGGAUUCUGAUUCCGGAGGUCGAGGUUCGAAUCCUCGUACCCCAGCCA
               4  3   5       3  3            4                          530  43 G
              45 34           5 34           44                         5303  33 G
              54 43           4 44          443                         3033  33 G
               4 3    5       3 4          443                          0334  3  G
                3    5    1343 4          4                          4 0    43   U
                         0                 3        4                 0          A
                   55  50 3243 3    5     3      5 4   5      5      0 5         U
                     3 03   34    55     3      5 44  4      5       1           C
                      0  0    1 55     53           23    5 5    5443 4          G
                     0 40   51     5              52     5   5       4           C
                       04   1     5               25    5    5       5           C
                         3           45                    5    4  5             A
                         3           5                    54    5  4             A
                         3    2 55                        4 5     4 5 5       5  G
                            42    3            5        54   3       5       5   C
                              4 43            5      5    5 3    45      5 5   4 G
                              2 35           5       4    535   45      5 5   4  G
                             3 3                   5  5545     4       5     4   U
                                                    4                 5          A
                                        5           3              4  5          A
                                44     5      4     10     44     4 5     54   1 G
                                 5    5      4      03    445    4 5  4  54   1  G
                                                5 10     4           4       1   C
                                           5        3     4           5          A
                                               40 33    54           5           C
                                         4     03 3   4 4    5                   C
                                        4    40      4      5     2        4   5 G
                                       4    305     43     5     2        4   5  G
                                        5  20       3     5     2                A
                                         32      533     5     2                 U
                                          1      345    5     2                  U
                                               53 5          2                   C
                                               3 5    5      4                   U
                                              3      5     33                    G
                                                    5     33                     A
                                                   5    43                       U
                                                  5    433                       U
                                                   5  4 34                       C
                                                  55  3 4                        C
                                                     3            4 5     52   5 G
                                                    34           45 1    523  5  G
                                                                4  1  4          A
                                                                53155 4 2134   4 G
                                                                31455 2 134   4  G
                                                              13     5 1     4   U
                                                             1       0           C
                                                           51    4450 5   54     G
                                                           1    4  0             A
                                                            4    3045    3       G
                                                                30455   3        G
                                                              31       3         U
                                                               1       5         U
                                                                     3           C
                                                                 55 3            G
                                                                   3             A
                                                                                 A
                                                                                 U
                                                                                 C
                                                                                 C
                                                                                 U
                                                                                 C
                                                                                 G
                                                                                 U
                                                                                 A
                                                                                 C
                                                                                 C
                                                                                 C
                                                                                 C
                                                                                 A
                                                                                 G
                                                                                 C
                                                                                 C
                                                                                 A

        green basepairs are kept as decision variables.
```
This triangular matrix gives you the probability that two nucleotides are canonically paired. For example, a 5 means there is a 10^-5 probability that the two corresponding bases are paired. For this reason, Biorseo does not consider all possible basepairs, but only the most probable ones, highlighted in green. You can set the threshold by using the `--theta` parameter in the command line options, default is 0.001.
*Note : these probabilities are computed using ViennaRNA's [vrna_pfl_fold](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__part__func__window.html) routine, which considers windows of size 100 and only allows basepairs within these windows of size 100. #TODO: allow user to customize window size.*
### PREPARING ILP DECISION VARIABLES AND CONSTRAINTS

Then the program defines the variables : 
```
Defining problem decision variables...
 > Legal basepairs : 0-67 1-66 2-65 3-64 4-18 4-63 5-17 5-62 6-16 6-19 6-61 7-15 7-61 8-14 8-17 8-22 8-44 9-13 9-16 9-21 9-43 10-15 10-20 10-42 13-22 14-21 16-22 20-44 20-45 20-71 21-44 21-70 22-42 22-43 22-69 24-40 25-39 26-38 26-58 27-37 27-57 28-35 28-36 28-56 29-34 29-55 30-34 30-54 31-53 39-67 40-60 40-66 41-59 42-58 42-64 42-65 43-57 43-62 43-64 44-54 44-63 45-53 45-61 46-52 46-60 47-51 47-59 48-58 49-57 50-55 51-55 
 > The possible stacks of two base pairs (i,j) and (i+1,j-1) : 0-67 1-66 2-65 3-64 4-18 4-63 5-17 5-62 6-16 7-15 8-14 8-17 8-22 8-44 9-16 9-21 9-43 13-22 20-45 20-71 21-44 21-70 24-40 25-39 26-38 26-58 27-37 27-57 28-35 28-56 29-55 30-54 39-67 40-60 41-59 42-58 42-65 43-64 44-54 45-53 45-61 46-52 46-60 47-59 48-58 
```
This summarizes the possible basepairs and stacks, depending on the previous probability matrix.

```
 > Looking for insertion sites...
 > Ignoring motif "1Z7F.A.2", hairpin (terminal) loops must be at least of size 3 !
 > Ignoring motif "3L26.C.2", hairpin (terminal) loops must be at least of size 3 !
 ...
 ```
 Depending on you module data source, you may then get messages about your modules, informing you if they could be placed in the sequence or not, sometimes because of errors. The results are then summarized : 

 ```
 > 262 candidate motifs on 4695 (53 ignored motifs), 
          50 insertion sites kept after applying probability threshold of 0.001
        > Allowed candidate insertion sites:
                > 3CUL.D.6      3CUL.D.6 ( 48-58 )         1 components: ........... basepairs:
                > 3IZF.A.2      3IZF.A.2 ( 24-26 38-40 )   2 components: ... ... basepairs:
                ...
```
*Note : for modules that are loops (DESC or BGSU's ones), the loop closing basepairs ARE considered even if they do not display here.*

This means, we had 4695 modules in the dataset, 53 of them were invalid, and we could potentially place 262 of them in the input sequence (not at the same time !). The optimisation program will now decide which of them we keep.

The program then defines constraints for the integer-linear-program.
```
        > 71 + 45 + 93 (yuv + xuv + Cpxi) decision variables are used.
        > ensuring there are at most 1 pairing by nucleotide...
            ...
        > ensuring that the stacks are correct...
            ...
        > forbidding lonely basepairs...
            ...
        > forbidding basepairs inside included motif's components...
            ...
        > forbidding component overlap...
            ...
        > ensuring that motives cannot be partially included...
            ...
        > forcing basepairs imposed by a module insertion...
            ...
    A total of 492 constraints are used.
```

### SOLVING

The program starts by solving the two objectives independantly : 
```
Solving...
        > Solution status: objective values (3.349, 15.0427)
        > Solution status: objective values (0, 19.1801)
Best solution according to objective 1 :(((((((.((...)).))...((.((.((.......)).))..))((.((.......)).)).)))))..... + 1VOY.B.63 + 1XJR.A.3 + 1YI2.0.158 + 2LA5.A.1 + 2RD2.B.3     3.3489961  15.0426987
Best solution according to objective 2 :((((((((((...)))....(((.((((([[[....)))))..)))((((...]]].))))))))))).....       0.0000000       19.1801369
```
Then, it will limit iteratively the criterion space (approach similar to the epsilon-constraint method) using constraints, and find new solutions to add to the Pareto front, or sometimes discard some dominated ones. This can last several hundreds or thousand iterations.

```
...

Solving objective function 2, on top of 1.027879251: Obj1  being in [1.027889251, 3.34900609441]...
        > Solution status: objective values (1.48611998869, 18.2887755166), not dominated.
        > adding structure to Pareto set :      ((((.(((((...)))....(((.((((([[[....)))))..)))((((...]]].)))))).))))..... + 1VOY.B.61 + 1XJR.A.3        1.4861200       18.2887755

...
```
### RESULTS

Finally, you get the Pareto set:
```
---------------------------------------------------------------
Whole Pareto Set:

        GGGGUAUCGCCAAGCGGUAAGGCACCGGAUUCUGAUUCCGGAGGUCGAGGUUCGAAUCCUCGUACCCCAGCCA
        ((((((((((...)))....(((.((((([[[....)))))..)))((((...]]].)))))))))))..... + 1XJR.A.3    0.7124144       19.1801369
                              ---               ----                                    1XJR.A.3 ( 22-24 40-43 )


        GGGGUAUCGCCAAGCGGUAAGGCACCGGAUUCUGAUUCCGGAGGUCGAGGUUCGAAUCCUCGUACCCCAGCCA
        ((((((((((...)))....(((.(((((.......)))))..)))((((.......)))))))))))..... + 1XJR.A.3 + 2RD2.B.3 1.0278793       19.1770336
                              ---               ----                                    1XJR.A.3 ( 22-24 40-43 )
                                    ---------                                           2RD2.B.3 ( 28-36 )


        GGGGUAUCGCCAAGCGGUAAGGCACCGGAUUCUGAUUCCGGAGGUCGAGGUUCGAAUCCUCGUACCCCAGCCA
        ((((.(((((...)))....(((.((((([[[....)))))..)))((((...]]].)))))).))))..... + 1VOY.B.61 + 1XJR.A.3        1.4861200       18.2887755
           ---                                                        ---               1VOY.B.61 ( 3-5 62-64 )
                              ---               ----                                    1XJR.A.3 ( 22-24 40-43 )


        GGGGUAUCGCCAAGCGGUAAGGCACCGGAUUCUGAUUCCGGAGGUCGAGGUUCGAAUCCUCGUACCCCAGCCA
        ((((.(((((...)))....(((.(((((.......)))))..)))((((.......)))))).))))..... + 1VOY.B.61 + 1XJR.A.3 + 2RD2.B.3     1.8015849       18.2856722
           ---                                                        ---               1VOY.B.61 ( 3-5 62-64 )
                              ---               ----                                    1XJR.A.3 ( 22-24 40-43 )
                                    ---------                                           2RD2.B.3 ( 28-36 )


        GGGGUAUCGCCAAGCGGUAAGGCACCGGAUUCUGAUUCCGGAGGUCGAGGUUCGAAUCCUCGUACCCCAGCCA
        ((((.(((((...)))....(((.(((((.......)))))..)))(((.........))))).))))..... + 1VOY.B.61 + 1XJR.A.3 + 2RD2.B.3 + 3CUL.D.6  2.0906497       17.3138502
           ---                                                        ---               1VOY.B.61 ( 3-5 62-64 )
                              ---               ----                                    1XJR.A.3 ( 22-24 40-43 )
                                    ---------                                           2RD2.B.3 ( 28-36 )
                                                        -----------                     3CUL.D.6 ( 48-58 )


        GGGGUAUCGCCAAGCGGUAAGGCACCGGAUUCUGAUUCCGGAGGUCGAGGUUCGAAUCCUCGUACCCCAGCCA
        ((((.(((((...)))....(((.((.(([[[....)).))..)))((((...]]].)))))).))))..... + 1FFK.0.86 + 1VOY.B.61 + 1XJR.A.3    2.2598256       17.2909202
                                 ---         ---                                        1FFK.0.86 ( 25-27 37-39 )
           ---                                                        ---               1VOY.B.61 ( 3-5 62-64 )
                              ---               ----                                    1XJR.A.3 ( 22-24 40-43 )


and several more...
```
There can be up to tens of solutions in here, depending mostly on your sequence length and average module size. For each solution, you get a secondary structure, and the corresponding insertion sites of some modules that match the sequence (they may have several components, which means several sequence portions are concerned by the same module). You get the insertion coordinates for each module. Also note the two criterion scores of the solution next to the secondary structure and the list of its modules.

The solutions are summarized in the result file, without the insertion sites locations : 

```
$ cat data/MyOutput.biorseo

PDB_00376
GGGGUAUCGCCAAGCGGUAAGGCACCGGAUUCUGAUUCCGGAGGUCGAGGUUCGAAUCCUCGUACCCCAGCCA
((((((((((...)))....(((.((((([[[....)))))..)))((((...]]].)))))))))))..... + 1XJR.A.3    0.7124144       19.1801369
((((((((((...)))....(((.(((((.......)))))..)))((((.......)))))))))))..... + 1XJR.A.3 + 2RD2.B.3 1.0278793       19.1770336
((((.(((((...)))....(((.((((([[[....)))))..)))((((...]]].)))))).))))..... + 1VOY.B.61 + 1XJR.A.3        1.4861200       18.2887755
((((.(((((...)))....(((.(((((.......)))))..)))((((.......)))))).))))..... + 1VOY.B.61 + 1XJR.A.3 + 2RD2.B.3     1.8015849       18.2856722
((((.(((((...)))....(((.(((((.......)))))..)))(((.........))))).))))..... + 1VOY.B.61 + 1XJR.A.3 + 2RD2.B.3 + 3CUL.D.6  2.0906497       17.3138502
((((.(((((...)))....(((.((.(([[[....)).))..)))((((...]]].)))))).))))..... + 1FFK.0.86 + 1VOY.B.61 + 1XJR.A.3    2.2598256       17.2909202
((((.(((((...)))....(((.((.((.......)).))..)))((((.......)))))).))))..... + 1VOY.B.61 + 1XJR.A.3 + 2LA5.A.1 + 2RD2.B.3  2.5752905       17.2878169
((((.(((((...)))....(((.((.((.......)).))..)))(((.........))))).))))..... + 1VOY.B.61 + 1XJR.A.3 + 2LA5.A.1 + 2RD2.B.3 + 3CUL.D.6       2.8643553       16.3159949
(((.((((((...))))....((.((.(([[[....)).))..))((.((...]]].)).)))).)))..... + 1FFK.0.76 + 1XJR.A.3 + 2JLT.B.1 + 2LA5.A.1  3.0335312       15.1786201
(((.((((((...))))....((.((.((.......)).))..))((.((.......)).)))).)))..... + 1XJR.A.3 + 2JLT.B.1 + 2LA5.A.1 + 2RD2.B.3 + 430D.A.1        3.3489961       15.1755168
(((.((((((...))))....((.((.(([[[....)).))..))((.((...]]].)).)))).)))..... + 1N36.A.99 + 1XJR.A.3 + 1YI2.0.158 + 2LA5.A.1        3.0335312       15.1786201
((((.(((((...)))....(((.((.((.......)).))..)))(((.........))))).))))..... + 1SER.T.3 + 1VOY.B.61 + 1XJR.A.3 + 2LA5.A.1 + 2RD2.B.3       2.8643553       16.3159949
((((.(((((...)))....(((.((.((.......)).))..)))((((.......)))))).))))..... + 1FFK.0.86 + 1VOY.B.61 + 1XJR.A.3 + 2RD2.B.3 2.5752905       17.2878169
((((.(((((...)))....(((.((.(([[[....)).))..)))((((...]]].)))))).))))..... + 1VOY.B.61 + 1XJR.A.3 + 2LA5.A.1     2.2598256       17.2909202
((((.(((((...)))....(((.(((((.......)))))..)))(((.........))))).))))..... + 1SER.T.3 + 1VOY.B.61 + 1XJR.A.3 + 2RD2.B.3  2.0906497       17.3138502
```
You still get the secondary structure, the name of the modules, and the two scores.
*Note : on linux/mac, if you ran Biorseo using Docker, the produced file MyResult.biorseo is owned by root. You may need to change permissions or ownership to access it with a regular user account.*

