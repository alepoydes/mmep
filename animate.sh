#!/bin/bash
#rm -f temp/*
#octave octave/animate.m $1
convert -delay 50 -loop 0 temp/*png $1.gif
#ffmpeg -i "./temp/%05d.png" -y $1.mpeg
#mencoder mf://output/*.png -mf w=640:h=480:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $1.avi
#rm -f temp/*