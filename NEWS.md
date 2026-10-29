Unreleased
----------

Notable changes:

 * Bugfix: in `--meth` mode indels were not fully left aligned. A reference C
   aligned to a read T scores a match under the methylation matrix, so shifting
   an indel across such a pair costs nothing, but the left alignment in
   `mb_fix_cigar` compared bases for equality and stopped there. Deletions
   therefore came to rest to the right of the position SAM and VCF convention
   asks for. The shift test now treats two bases as interchangeable when both
   match the base they swap against, which is equality outside `--meth`, so
   nothing changes there. On a 2 Mb test window this moved 3,891 indels left and
   none right, left the alignment score and `NM` of every read unchanged, and
   raised exact CIGAR agreement with bwa-meth (which uses bwa-mem for alignment) on indel-bearing reads from 75.9% to 88.8% with no reads newly disagreeing.

Release 0.7-r421 (6 August, 2026)
---------------------------------

Notable changes:

 * Bugfix: misassgined "properly paired" SAM flag (#60). After this fix, flag
   0x2 is only applied if the primary hits land on the same chromosome. The
   flag for a supplementary hit follows the primary hit for consistency with
   bwa-mem.

 * Bugfix: unrecognized command-line options may be reported at wrong
   positions.

(0.7: 6 August 2027, r421)



Release 0.6-r416 (30 JUly, 2026)
--------------------------------

Notable changes:

 * Bugfix: option `-I` was not working (#57 and #58)

(0.6: 30 July 2026, r416)



Release 0.5-r414 (25 July 2026)
-------------------------------

Notable changes:

 * Improvement: theoretically better mapping quality (mapQ) for paired-end
   reads. Bwa-mem overestimates mapQ in a rare case, which is an oversight.
   The new release fixes the issue. On real-world SNP calling, the change
   increases precision but reduces sensivity. See 98e89de for details.

 * New feature: option `--mmap` for loading index as memory-mapped files (#55).
   It reduces index loading time but may slow down mapping depending on the
   system.

 * New feature: option `-I` to specify the insert size distribution without
   inferring from data (#57).

 * New feature: option `--outs` to skip poor hits in output or in the XA tag
   (#56).

 * Bugfix: APIs did not support BS-seq (#53).

(0.5: 25 July 2026, r414)



Release 0.4-r400 (12 July 2026)
-------------------------------

Notable changes:

 * Improvement: penalize heavily clipped alignment when rescoring. Older
   minibwa may produce a clipped primary alignment together with an end-to-end
   supplementary alignment with more mismatches/gaps. While often neither
   alignment is correct, having them both is confusing. This happens very
   rarely.

 * Bugfix: subcommand name was missing from the PG-CL tag in output SAM (#46).
   Also fixed a few typos in command-line help messages (#48 and #49).

(0.4: 12 July 2026, r400)



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
