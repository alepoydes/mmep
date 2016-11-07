#!/bin/sh
mkdir -p tmp/
bin/mepd $@ 
gnuplot tmp/mep.gnuplot
gnuplot tmp/energy.gnuplot