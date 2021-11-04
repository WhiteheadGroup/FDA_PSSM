# Code to generate FR Scores

See published article at https://www.frontiersin.org/articles/10.3389/fimmu.2021.728694/full

Finds and scores all framework mutations from input antibody file (csv format). Outputs normalized FR scores as a scatterplot. Verbose mode prints full pairwise alignment of each antibody. The output_mutations option creates a csv with all individual mutation scores. The output_csv option creates a csv with overall FR scores for each analyzed antibody.


- Input sequences are entered into FDA_Abs.csv

- Normalization constants are entered into normalization.csv

- Scatterplot output as svg file named FDA_Abs.svg

- Results from a sample run are included
