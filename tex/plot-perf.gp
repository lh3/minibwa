set t po eps co so "Helvetica,20"

set out "fig-perf.eps"
set size 1.9,1.6
set multiplot layout 2,3

set style data histogram
set style histogram cluster gap 2
set style fill solid border -1
set boxwidth 0.85
set xtic rotate by -45 scale 0
set border 3
set bmargin 5.2

# Top row: speed
set origin 0,0.8
set size 1,0.8
set yran [0:160]
set ylab "Speed (Gbp/hr over 32 threads)"

# Short reads
set key top left
plot "<cat perf.tsv|grep WGS" using 3:xtic(2) ti "WGS" ls 2, \
	"<cat perf.tsv|grep SBX" using 3 ti "SBX" ls 3, \
	"<cat perf.tsv|grep Hi-C" using 3 ti "Hi-C" ls 4

# Long reads
set origin 0.95,0.8
set size 0.6,0.8
unset ylab
set key top right
plot "<cat perf.tsv|grep HiFi" using 3:xtic(2) ti "HiFi" ls 1, \
	"<cat perf.tsv|grep ONT" using 3 ti "ONT" ls 5 lc "#d95f0e"

# Meth reads
set style histogram cluster gap 1
set origin 1.5,0.8
set size 0.4,0.8
plot "<cat perf.tsv|grep BS-seq" using 3:xtic(2) ti "BS-seq" ls 6

# Bottom row: memory
set yran [0:140]
set ylab "Peak memory (GB)"

# Short reads
set style histogram cluster gap 2
set size 1,0.8
set key top left
plot "<cat perf.tsv|grep WGS" using 4:xtic(2) ti "WGS" ls 2, \
	"<cat perf.tsv|grep SBX" using 4 ti "SBX" ls 3, \
	"<cat perf.tsv|grep Hi-C" using 4 ti "Hi-C" ls 4

# Long reads
unset ylab
set key top right
set origin 0.95,0
set size 0.6,0.8
plot "<cat perf.tsv|grep HiFi" using 4:xtic(2) ti "HiFi" ls 1, \
	"<cat perf.tsv|grep ONT" using 4 ti "ONT" ls 5 lc "#d95f0e"

# Meth reads
set style histogram cluster gap 1
set origin 1.5,0
set size 0.4,0.8
plot "<cat perf.tsv|grep BS-seq" using 4:xtic(2) ti "BS-seq" ls 6

unset multiplot
