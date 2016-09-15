#!/bin/sh
bin/minf $@ -p | gnuplot -p 2>/dev/null
gnuplot tmp/state.gnuplot