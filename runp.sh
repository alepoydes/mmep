#!/bin/bash
type sem >/dev/null 2>&1 || { echo -e >&2 "Install required program:\nsudo apt-get install parallel"; exit 1; }



PREFIX="tmp/"
PARAM=$1
shift 1

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
	sem -j +0 --id $$ $@ -D "${path}"
	gnuplot "${path}/"*.gnuplot
done < ${PARAM}
sem --wait --id $$