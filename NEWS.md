Unreleased
----------

Notable changes:

 * New feature: BAM input. Minibwa reads BAM wherever it reads FASTA/FASTQ,
   detecting the format from the file contents. Both unaligned and aligned BAM
   work; for aligned input the old alignments are discarded and secondary and
   supplementary records are skipped. Paired-end mode is turned on
   automatically when the first record is flagged as paired.

 * Bugfix: in `--meth` mode the bisulfite conversion was chosen from a read's
   position within its fragment, so an unpaired read 2 -- which a BAM can carry
   but a pair of FASTQ files cannot express -- was converted as if it were
   read 1 and usually failed to map. With a BAM input the SAM FLAG now decides,
   and the segment a read belongs to is carried through to the output FLAG, so
   an unpaired read 2 is still reported as read 2 rather than as single-end.

 * New feature: paired-end input is now validated. A coordinate-sorted BAM is
   rejected in paired-end mode instead of being mapped as single-end, and two
   query files whose reads are out of order are reported rather than silently
   mispaired.

 * New feature: with a BAM input, `@RG`, `@PG` and `@CO` header records and
   per-record auxiliary tags such as `RX` are carried over to the output.
   `@HD` and `@SQ` are replaced, as contigs come from the index, and tags
   describing an alignment are dropped because minibwa recomputes them. Use
   `--bam-tags=no` to copy no tags.



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
