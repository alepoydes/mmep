#!/bin/sh
bin/mind $@ -p | gnuplot -p 2>/dev/null
gnuplot tmp/state.gnuplot