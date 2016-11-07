#!/bin/sh
basedir=$1
mkdir -p ${basedir}
shift 1
output=${basedir}/rate.txt
echo "Transition rates" > "${output}"
for d in ${basedir}/*; do
	echo "===========================================" >> "${output}"
	echo "Processing $d/mep.oct" >> "${output}"
	echo "===========================================" >> "${output}"
	./rate.py "$d/mep.oct" $@ >> "${output}" 2>> "${output}"
done