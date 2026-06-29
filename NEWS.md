Release 0.3-r391 (28 June 2026)
-------------------------------

Notable changes:

 * New feature: added the `mem` subcommand to mimic the bwa-mem command-line
   interface (CLI). Most input/output options and commonly used options are
   retained; unsupported or incompatible bwa-mem options are silently ignored.
   For commonly used bwa-mem command lines, replacing "bwa mem" with "minibwa
   mem" often gives comparable outputs.

 * New feature: added option `-H` to inject arbitrary header lines.

 * New feature: added the `XA` tag to encode secondary alignments if there are
   not too many (#25).

 * Bugfix: SAM header was outputted to stdout even if an output file is
   specified.

 * Bugfix: seed sorting was incorrectly implemented (#39).

 * Bugfix: compiled mimalloc missed the `aligned_alloc` function in C11 (#30).

 * Bugfix: methylation index was incorrectly generated under the low-memory
   indexing mode (#42).

(0.3: 28 June 2026, r391)



Release 0.2-r370 (18 June 2026)
-------------------------------

This release fixed several minor bugs:

 * Bugfix: FASTA comments not parsed correctly during indexing (#22). It is
   recommended to reindex the reference genome with the new release.

 * Bugfix: double TABs ahead of each MD tag (#27)

 * Bugfix: a small memory leak during index loading (#26)

 * Bugfix: compilation error on ARM Linux (#24)

(0.2: 18 June 2026, r370)



Release 0.1-r363 (15 June 2026)
-------------------------------

First public release.

(0.1: 15 June 2026, r363)
