## Overview

Minibwa combines bwa-mem seeding with minimap2 chaining and base alignment. It
aims to replace bwa-mem for short-read alignment with better performance at
comparable accuracy. It also works for accurate long read alignment to some
extent. The following table summarizes major components in minibwa and their
relationship with bwa-mem, ropebwt3 and minimap2:

|Component|Source  | Modifications |
|:--------|:-------|:--------------|
|FM-index |bwa-mem | Near identical |
|SMEM     |ropebwt3| Reimplemented and batched |
|SA query |ropebwt3| Batched differently |
|Seeding  |bwa-mem | Different algorithm |
|Chaining |minimap2| Adapted for variable-length seeds |
|Alignment|minimap2| Same ksw2 routines but different gluing code |
|Pairing  |bwa-mem | With pre-alignment filtering |

## Similarity to bwa-mem and minimap2

Similarity to bwa-mem:

 * Almost the same data structure for FM-index. Nonetheless, at the API level,
   minibwa uses half-closed-half-open suffix array (SA) intervals, similar to
   ropebwt3, but bwa-mem uses closed intervals. Minibwa also increases the SA
   sampling rate from 1/32 to 1/16.

Similarity to minimap2:

 * Minibwa chaining is adapted from minimap2, with modifications for
   variable-length seeds in minibwa. Bwa-mem tree-based chaining is no longer
   used.

 * Like minimap2, minibwa uses ksw2 for base alignment. Bwa-mem extension code
   is removed.

## Differences from bwa-mem and minimap2

 * Minibwa constructs FM-index with libsais by default. This library supports
   multi-threading (requiring OpenMP). It is much faster at the cost of more
   memory. The old bwa-mem BWT construction algorithm is still available.

 * Although the memory layout of FM-index is near identical, the index file
   formats are very different now. Minibwa introduces `l2bit` which is similar
   to UCSC's 2-bit format but supports contigs longer than 32-bit integers.

 * Minibwa uses Travis Gagie's algorithm for SMEM finding (function
   `mb_bwt_smem`). This algorithm was first implemented in ropebwt3. Although
   for standard SMEM finding in short reads, it is only slightly faster than
   the original algorithm, the new method is simpler and enables further
   optimization:

   - Minibwa caches the SA intervals of all 10-mers (`mb_bwt_cache`) which
     reduces random memory access. The old SMEM algorithm is incompatible with
     this pattern.

   - Minibwa batches multiple sequences for SMEM finding to hide latency with
     prefetch (`mb_bwt_smem_batch`), speeding up SMEM finding by 2.5 folds.
     Coding agents offerred great help in debugging the complex logic involved.
     It is possible to accelerate the old SMEM algorithm with latency hiding
     but this will be much more complex.

   - Minibwa uses two rounds of the SMEM algorithm for seeding
     (`mb_seed_intv`) and thus benefits from latency hiding throughout the
     seeding step (`mb_seed_intv_batch`). Bwa-mem does not use this strategy
     because with the old algorithm, the second round would be inefficient.

 * When retrieving SA values, minibwa also batches multiple BWT positions
   (`mb_bwt_sa_batch`). Bwa-mem2 uses the same trick.

 * Minibwa slightly speeds up SSE-based alignment originated from minimap2. The
   changes will be merged back to minimap2 later.

 * During mate rescue, minibwa uses a small 7-mer hash table to filter out poor
   alignment without full Smith-Waterman (`mb_ungap`). The algorithm is
   inspired by the Hough transform for line finding.

 * Minibwa reduces alignments in highly repetitive regions because short reads
   in these regions cannot be aligned well anyway. On simulated data, the
   accuracy remains similar.

 * Minibwa natively aligns directional BS-seq reads.

## BAM input and aux tags

Minibwa reads BAM directly (`bam.c`). BGZF is a multi-member gzip stream, so
zlib's `gzread` inflates it and only the binary record layout has to be decoded;
this keeps minibwa free of an htslib dependency, at the cost of single-threaded
decompression. Aligned BAM is accepted as well as unaligned: records flagged
0x10 are reverse complemented back to the original read, and secondary and
supplementary records are skipped so that each read is mapped once (a
hard-clipped supplementary record would otherwise contribute a truncated
sequence).

Header records other than `@HD` and `@SQ` are carried over, `@RG` deduplicated
by ID across inputs; contigs come from the index, so the input `@SQ` lines are
replaced. Minibwa's own `@PG` chains to the last input `@PG` via `PP:`, and its
`ID:` is uniquified to `minibwa.1`, `minibwa.2` and so on when the input already
contains one, since `@PG` identifiers must be unique.

In `--meth` mode the in-silico conversion follows the SAM FLAG when the input
is a BAM: read 1 gets C-to-T and read 2 G-to-A regardless of whether the mate is
present. FASTA/FASTQ carries no flag, so there the position within the fragment
is still used. This matters because mixed single- and paired-end input is
allowed, so an unpaired read 2 is legal and must not be converted as a read 1.
The same rule sets the output FLAG: a read the input marked as read 1 or read 2
keeps that bit even when its mate is absent, with 0x8 set to say the mate is not
there.

Per-record aux tags are carried over except those listed in `bam_aux_drop[]`
(`bam.c`), which are dropped for one of two reasons:

 1. Minibwa writes the tag itself. [SAM][samv1] requires that "within each
    alignment line, no `TAG` may appear more than once", so passing the input copy
    through would produce an invalid record. Minibwa writes `NM`, `AS`, `ms`
    and `md` on every aligned record, plus `MC`, `MQ`, `SA`, `XA`, `CG` and
    `n2`, and `cs`, `ds` or `MD` under `-b`. This group is belt and braces:
    `sam_write_aux()` appends the carried-over tags last and skips any name
    already present on the record, so validity does not depend on the list
    being complete. A gap costs one tag, not a malformed record.

 2. The tag describes an alignment made by some other program and says nothing
    true about the alignment minibwa is about to produce.

`RG` is deliberately not in the list: it is a property of the read rather than
of the alignment, so it is carried over. `-R` overrides it: an explicit read
group on the command line replaces whatever the records carried, and the input
`@RG` header lines are dropped because nothing refers to them any more.

The complete list of dropped tags and the logic behind tag handling is in [bam.c](bam.c)
(`bam_aux_drop[]`) with the writing logic in [format.c](format.c)
(`sam_write_aux()`), and minibwa's own tag definitions are in [minibwa.1](minibwa.1) "OUTPUT FORMATS"
and in the code ([cs.c](cs.c) `mb_write_cs_ds`, [cs.c](cs.c) `mb_write_MD`, [format.c](format.c) `n2`).

[samv1]: https://samtools.github.io/hts-specs/SAMv1.pdf
