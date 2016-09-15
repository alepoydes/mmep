#!/bin/sh
bin/mepd $@ 
gnuplot tmp/mep.gnuplot
gnuplot tmp/energy.gnuplot