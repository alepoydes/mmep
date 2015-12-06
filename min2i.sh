#!/bin/sh
./min2q $@ -p | gnuplot -p 2>/dev/null
gnuplot fields/mep.gnuplot