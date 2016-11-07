#!/bin/bash
OUTPUT=tmp/batch.gnuplot
TMP=/tmp/$$
ENERGYFILE=tmp/energy.txt
DISTANCEFILE=tmp/distance.txt
OPTIMIZE="bin/mepf -o"
echo "# $@" > ${OUTPUT}
DESC=$1
SUFF=$2
shift 2

mkdir -p tmp/

echo -n "" > ${ENERGYFILE}
echo -n "" > ${DISTANCEFILE}
echo 'set title "Energy variation"' >> ${OUTPUT}
echo 'unset key' >> ${OUTPUT}
echo 'load "moreland.pal"' >> ${OUTPUT}
echo 'set cblabel "energy"' >> ${OUTPUT}
echo 'set xlabel "along path"' >> ${OUTPUT}
echo 'set ylabel "parameter"' >> ${OUTPUT}
echo 'set terminal pngcairo' >> ${OUTPUT}
echo 'set output "tmp/batch.png"' >> ${OUTPUT}

echo '$data << EOD' >> ${OUTPUT}
let "a=1"
while IFS='' read -r line || [[ -n "$line" ]]; do
	eval "${line}"
	#envsubst < ${DESC} 
	envsubst < ${DESC} | ${OPTIMIZE} $@ > ${TMP}
	cat ${TMP} >> ${OUTPUT}
	sed -n '1p' ${TMP} >> ${ENERGYFILE}
	sed -n '2p' ${TMP} >> ${DISTANCEFILE}
	gnuplot tmp/mep.gnuplot
	gnuplot tmp/energy.gnuplot
	mv -f tmp/mep.gnuplot tmp/mep.${a}.gnuplot
	mv -f tmp/energy.gnuplot tmp/energy.${a}.gnuplot
	mv -f tmp/mep.gif tmp/mep.${a}.gif
	mv -f tmp/mep.oct tmp/mep.${a}.oct
	mv -f tmp/energy.png tmp/energy.${a}.png
	let "a=a+1"
done < ${SUFF}
echo "EOD" >> ${OUTPUT}

echo 'set view map' >> ${OUTPUT}
echo 'splot "$data" matrix with image' >> ${OUTPUT}

gnuplot5 tmp/batch.gnuplot

rm -f ${TMP}