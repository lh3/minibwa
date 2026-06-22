[![Bioconda](https://img.shields.io/conda/dn/bioconda/minibwa.svg?style=flag&label=Bioconda)](https://anaconda.org/bioconda/minibwa)
[![Homebrew](https://img.shields.io/homebrew/v/minibwa)](https://formulae.brew.sh/formula/minibwa)
[![preprint](https://img.shields.io/badge/arXiv-2606.15357-blue)](https://arxiv.org/abs/2606.15357)

## Getting Started
```sh
git clone https://github.com/lh3/minibwa
cd minibwa && make

# with test data
./minibwa index test/chrM-human.fa.gz chrM-human              # index the genome
./minibwa map chrM-human test/chrM-read_?.fa.gz > aln.sam     # align and output in SAM

# other examples without test data
minibwa map -ft16 ref.index long-read.fq > aln.paf            # align long reads
minibwa map --hic ref.index reads.interleaved.fq > aln.sam    # align Hi-C short reads

# align *directional* bisulfite sequencing (BS-seq) reads
minibwa index --meth -t8 ref.fa                               # generate BS-seq index
minibwa map --meth ref.fa read1.fq read2.fq > aln.sam         # map BS-seq reads
```

## Introduction

Minibwa aligns short reads against a reference genome. It is the successor of
[bwa-mem][bwa] with a different algorithm. Minibwa is over three times as fast as the
original bwa-mem and twice as fast as [bwa-mem2][bwa-mem2] at comparable accuracy. While
minibwa works with accurate long reads, [minimap2][mm2] is more robust under high
error rate.

Minibwa is a hybrid of bwa-mem and minimap2: it indexes the genome with
Burrow-Wheeler Transform (BWT), finds variable-length seeds like bwa-mem, and
performs chaining and SIMD-based nucleotide alignment with the minimap2
algorithm. Minibwa speeds up bwa-mem2 further with additional prefetch for
seeding, new heuristics to skip unnecessary mate rescue and reduced effort in
highly repetitive regions where reads would often be wrongly mapped due to
structural changes anyway.

## Users' Guide

### Intended use cases

Minibwa is designed for mapping short reads and accurate long reads. It does
not support spliced alignment and has not been tuned for aligning long contigs.
For now, minibwa does not properly work with alternate contigs in the reference
genome. Please use a version of the reference without such contigs.

### Installation

Minibwa requires either SSE4.2 on x86 CPUs or NEON on ARM. It depends on
[zlib][zlib] installed on your system and also includes slightly modified
source code of [mimalloc][mimalloc] and [libsais][libsais] which optionally
uses OpenMP for multi-threading. You can build minibwa with
```sh
make             # automatically detect OpenMP and arm64 vs. x86_64
make omp=0       # disable multi-threading in libsais (no effect on mapping)
make gpl=0       # disable GPL'd code for low-memory BWT building (no effect on mapping)
make mimalloc=0  # disable mimalloc and use the system malloc+kalloc instead
```
This produces a single binary `minibwa` which you can copy to your `PATH`.

### Usage

Like bwa-mem, minibwa requires to index the genome before read alignment.

#### Indexing

You can index the reference genome with
```sh
minibwa index -t8 ref.fa     # index with 8 threads, using 18N RAM (N is the genome size)
minibwa index ref.fa prefix  # use a different index prefix instead of ref.fa
minibwa index -l ref.fa      # use less memory at the cost of performance
minibwa index --meth ref.fa  # generate BS-seq index
```
Minibwa generates two files: `ref.fa.l2b` for 2-bit encoded reference genome
sequences and `ref.fa.mbw` for BWT and sampled suffix array. In the `--meth`
mode, minibwa additionally generates `ref.fa.meth.mbw` for the BWT of the
3-base genome.

#### Mapping

By default, minibwa dynamically changes multiple internal parameters based on
individual read lengths. It works for both short and accurate long reads.
```sh
minibwa map -t8 ref.fa read1.fq read2.fq    # map paired-end reads and output SAM
minibwa map -ft8 ref.fa read.fa.gz          # map single-end or long reads; output PAF
minibwa map --hic ref.fa hic1.fq hic2.fq    # map Hi-C short reads
minibwa map --meth ref.fa read1.fq read2.fq # map BS-seq reads; requiring "index --meth"
```
Note in the default adaptive mode, `-g`/`-w`/`-W`/`-N`/`-m`/`-s` only changes
the short-read setting; the long-read setting is fixed. This mode is disabled
with `--adap=no` or when `-x sr` or `-x lr` is specified.

## Developers' Guide

Minibwa provides basic APIs for loading index and aligning reads.
[api-test/ex-one.c](api-test/ex-one.c) shows an example to align each read
independently; [api-test/ex-batch.c](api-test/ex-batch.c) aligns multiple reads
in batch, which is faster and also supports paired-end mapping.
[dev.md](dev.md) explains how minibwa differs from BWA-MEM and minimap2.

## Limitations

* Minibwa does not work with noisy long reads or spliced RNA-seq reads.
* Minibwa does not support undirectional bisulfite sequencing data.

[zlib]: https://zlib.net/
[mimalloc]: https://github.com/microsoft/mimalloc
[libsais]: https://github.com/IlyaGrebnov/libsais
[bwa]: https://github.com/lh3/bwa
[mm2]: https://github.com/lh3/minimap2
[bwa-mem2]: https://github.com/bwa-mem2/bwa-mem2
