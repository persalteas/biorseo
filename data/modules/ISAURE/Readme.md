The motif library used with --contacts is particular. It was provided by Isaure Chauvot de Beauchêne from the LORIA 
laboratory. These motifs are made up of RNA fragments linked to proteins.
==================================================================================================================

Several versions of these designs have been provided, but the most complete is the latest:'motifs_06-06-2021.json'
The current scripts were created based on this file, and doesn't work with the other older libraries.

There is also 2 benchmarks files also in json format : 'benchmark_16-06-2021.json' and 'benchmark_16-07-2021.json'.
It contains complete RNA sequences that bind to a protein, the first one contains only 33 RNA, and the second one 
contains 130 RNA.

The benchmark.dbn and benchmark.txt were created based on the 'benchmark_16-07-2021.json'. 
They are mostly used for the Isaure_benchmark.py script and scripts from the 'scripts' directory.

The motifs_final.json it obtains after executing the count_pattern.cpp script in 'script' directory on
the 'motifs_06-06-2021.json' motifs file.
This script count the number of "occurrences" of the motif. So we consider that if the sequence of motif A 
is included in motif B, then for each inclusion of B we also have an inclusion of A. And vice versa.

The motif library used by BiORSEO is the one in the 'bibliotheque_a_lire' directory. There should only be
the json file we wish to be used by BiORSEO for it's prediction. That's why you shouldn't put other type of file!






