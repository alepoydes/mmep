#!/bin/sh
bin/mepd $@ 
gnuplot fields/mep.gnuplot
gnuplot fields/energy.gnuplot