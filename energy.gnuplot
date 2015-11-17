set terminal png
set output '../fields/energy.png'
set ytics nomirror
set y2tics nomirror
set log y2
plot '../fields/energy.txt' using 1:2 with lines axes x1y1 title 'energy', '' using 1:3 with lines axes x1y2 title 'grad.', '' using 1:4 with lines axes x1y2 title 'orth. grad.'
#distance energy gradient