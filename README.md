#spades_stats

Calculating statistics for SPAdes assemblies based on the
coverage and length found in the fasta description line.

By Karin Lagesen | @karinlag 

## How to run

###Requirements

- Python 2.7
- Biopython

###Command line

```
usage: spades_stats.py [-h] [-d DIRECTORY] [-p STRING]

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Directory containing directories with scaffolds.fasta
                        in them
  -p STRING, --prefix STRING
                        Prefix string for output files

```

###Input

The input for this code is a directory where there are
spades assemblies. The code locates all files named
"scaffolds.fasta" under that directory and calculates 
statistics for that file.

###Output

Statistics is calculated for each file, and the results
are written to one output file. 

The output is output one line per scaffolds file, with
the following columns (first line is header line):

|Column           | Description 
|-------------------|------------ 
|Scaffolds_filename | name of input file, incl relative path
|N50              | N50 value for assembly
|\#contigs>=N50   | the fewest number of contigs whose sum makes up N50
|coverage_contigs | average coverage over the #contigs>=N50


Next comes five number statistics for length and coverage.
The five stats calculated are(in order):

- min: minimum value (length or coverage)
- max: maximum value (length or coverage)
- avg: average value (length or coverage)
- median: median value (length or coverage)
- std: standard deviation (length or coverage)

##Issues

Please report problems here: https://github.com/karinlag/spades_stats/issues

##License

Please read the LICENSE included in the repository

