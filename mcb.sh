#!/bin/bash
LOG="fields/mcb.log"
MC="bin/mcq"
date > ${LOG}
echo "mcb $@" >> ${LOG}
DESC=$1
SUFF=$2
shift 2

while IFS='' read -r line || [[ -n "$line" ]]; do
	eval "${line}"
	IFS=; components=(${line})
	echo $components
	for comp in ${components}; do
		IFS== read var val <<<"${comp}"
		export $var
	done
	filename="fields/${line}.bin"
	echo "${line} > ${filename}" >> ${LOG}
	#envsubst < ${DESC} 
	envsubst < ${DESC} | ${MC} $@ -o "${filename}"
done < ${SUFF}