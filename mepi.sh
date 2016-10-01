#!/bin/sh
bin/mepd $@ -p | gnuplot -p 2>/dev/null
gnuplot tmp/mep.gnuplot
gnuplot tmp/energy.gnuplot