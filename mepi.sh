#!/bin/sh
bin/mepq $@ -p | gnuplot -p 2>/dev/null
gnuplot fields/mep.gnuplot
gnuplot fields/energy.gnuplot