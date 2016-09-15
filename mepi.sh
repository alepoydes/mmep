#!/bin/sh
bin/mepf $@ -p | gnuplot -p 2>/dev/null
gnuplot tmp/mep.gnuplot
gnuplot tmp/energy.gnuplot