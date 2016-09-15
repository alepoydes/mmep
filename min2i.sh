#!/bin/sh
bin/min2q $@ -p | gnuplot -p 2>/dev/null
gnuplot tmp/mep.gnuplot