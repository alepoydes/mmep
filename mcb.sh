#!/bin/bash
LOG="tmp/mcb.log"
MC="bin/mcd"
date > ${LOG}
echo "mcb $@" >> ${LOG}
DESC=$1
SUFF=$2
shift 2

while IFS='' read -r line || [[ -n "$line" ]]; do
	IFS=';' read -ra components <<< "${line}" 
	filename="tmp/"
	for comp in "${components[@]}"; do
		IFS='=' read var val <<< "${comp}"
		export $var=$val
		echo "var" $var "val" $val "res" ${!var}
		filename="${filename}$var$val"
	done
	filename="${filename}.bin"
	echo "${line} > ${filename}" >> ${LOG}
	echo "run process"
	envsubst < ${DESC} | sem -j +0 --id $$ ${MC} $@ -o "${filename}"
done < ${SUFF}
sem --wait --id $$