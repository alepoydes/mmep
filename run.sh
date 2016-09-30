#!/bin/bash
PREFIX="tmp/"
PARAM=$1
shift 1

DESC=`echo $@ | tr -dc '[:alnum:]\n\r'`
PREFIX="${PREFIX}${DESC}/"
mkdir "${PREFIX}"

while IFS='' read -r line || [[ -n "$line" ]]; do
	IFS=';' read -ra components <<< "${line}" 
	path=$PREFIX
	for comp in "${components[@]}"; do
		IFS='=' read var val <<< "${comp}"
		export $var=$val
		path="${path}$var$val"
	done
	path="${path}/"
	echo "Populating " $path
	mkdir "${path}"
	$@ -D "${path}" || exit 1
	gnuplot "${path}/"*.gnuplot
done < ${PARAM}