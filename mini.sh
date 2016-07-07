#!/bin/sh
bin/mind $@ -p | gnuplot -p 2>/dev/null
gnuplot fields/state.gnuplot