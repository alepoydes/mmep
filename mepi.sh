#!/bin/sh
mkdir -p tmp/
bin/mepl $@ -p | gnuplot -p 2>/dev/null || exit 1
gnuplot tmp/mep.gnuplot 2>/dev/null
gnuplot tmp/energy.gnuplot 2>/dev/null
exit 0