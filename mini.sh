#!/bin/sh
mkdir -p tmp/
bin/minl $@ -p | gnuplot -p 2>/dev/null || exit 1
gnuplot tmp/state.gnuplot 2>/dev/null
exit 0