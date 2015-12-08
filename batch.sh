#!/bin/bash
OUTPUT=fields/batch.gnuplot
OPTIMIZE=./mepq
echo "# $@" > ${OUTPUT}
DESC=$1
SUFF=$2
shift 2

echo 'set title "Energy variation"' >> ${OUTPUT}
echo 'unset key' >> ${OUTPUT}
echo 'load "moreland.pal"' >> ${OUTPUT}
echo 'set cblabel "energy"' >> ${OUTPUT}
echo 'set xlabel "along path"' >> ${OUTPUT}
echo 'set ylabel "parameter"' >> ${OUTPUT}
echo 'set terminal pngcairo' >> ${OUTPUT}
echo 'set output "fields/batch.png"' >> ${OUTPUT}

echo '$data << EOD' >> ${OUTPUT}
let "a=1"
while IFS='' read -r line || [[ -n "$line" ]]; do
	eval "${line}"
	#envsubst < ${DESC} 
	envsubst < ${DESC} | ${OPTIMIZE} $@ >> ${OUTPUT}
	gnuplot fields/mep.gnuplot
	gnuplot fields/energy.gnuplot
	mv -f fields/mep.gnuplot fields/mep.${a}.gnuplot
	mv -f fields/energy.gnuplot fields/energy.${a}.gnuplot
	mv -f fields/mep.gif fields/mep.${a}.gif
	mv -f fields/energy.png fields/energy.${a}.png
	let "a=a+1"
done < ${SUFF}
echo "EOD" >> ${OUTPUT}

echo 'set view map' >> ${OUTPUT}
echo 'splot "$data" matrix with image' >> ${OUTPUT}

gnuplot5 fields/batch.gnuplot