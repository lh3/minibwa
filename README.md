## Getting Started
```sh
git clone https://github.com/lh3/minibwa
cd minibwa
make        # Or "make omp=0" if the compiler doesn't support OpenMP
./minibwa index test/chrM-human.fa.gz chrM-human              # index the genome
./minibwa map -a chrM-human test/chrM-read_?.fa.gz > aln.sam  # align and output in SAM
```

## Introduction

Minibwa aligns short and long reads against a reference genome. It is the
successor of bwa-mem for short-read alignment with a different algorithm.
Minibwa is 3-4 times as fast as the original bwa-mem and twice as fast as
bwa-mem2 at comparable accuracy. While minibwa works with accurate long reads,
it does not aim to replace minimap2. 

Minibwa is named after bwa-mem and minimap2 and is based on the source code of
both projects: it indexes the genome with Burrow-Wheeler Transform (BWT), finds
variable-length seeds as nested SuperMaximal Exact Matches (SMEMs) like
bwa-mem, and performs chaining and SIMD-based nucleotide alignment with the
minimap2 algorithm. Minibwa speeds up bwa-mem2 further with effective prefetch
for seeding, additional heuristics to skip unnecessary mate rescue and reduced
effort in highly repetitive regions where reads would often be wrongly mapped
anyway.
