#!/bin/bash
set -o verbose

type sem >/dev/null 2>&1 || { echo -e >&2 "Install required program:\nsudo apt-get install parallel"; exit 1; }

NUMJOBS=+0
#NUMJOBS=4

BASE="$(basename $1)_$(basename $2)_$(basename $3)"
DESC=$BASE

PREFIX="tmp/"
PARAM=$1
shift 1

PREFIX="${PREFIX}${DESC}/"
lock="${PREFIX}/.lock"

if [ -f ${lock} ]; then
	echo -e "\e[91mError:\e[0m Lock file ${lock} was not removed. Is directory in use?"
	exit 
fi

#shopt -s nullglob #safety needed so that globs return empty strings when no files are present
function finalize() {
	rm -f ${lock}
}

function killjobs() {
	rm -f ${lock}
    for s in ~/.parallel/semaphores/id-$$/*@*;
    do
       kill -15 -- -$(basename ${s%%@*})
    done
	exit 255
}

trap finalize EXIT
trap killjobs SIGINT

mkdir -p "${PREFIX}"
touch "${lock}"


dirs=()
while IFS='' read -r line || [[ -n "$line" ]]; do
	IFS=';' read -ra components <<< "${line}" 
	path="${PREFIX}"
	for comp in "${components[@]}"; do
		IFS='=' read var val <<< "${comp}"
		echo export "$var=$val"
		export "$var=$val"
		path="${path}$var$val"
	done
	path="${path}/"
	dirs+=("${path}")
	if [ "$(ls -A ${path} 2>/dev/null)" ]; then
		echo -e "\e[93mSkipping\e[0m non empty $path"
	else
		echo -e "\e[92mPopulating\e[0m ${path}"
		mkdir -p "${path}"
		sem -j ${NUMJOBS} --id $$ $@ -D "${path}"
	fi
done < ${PARAM}
sem --wait --id $$

for path in "${dirs[@]}"; do
	gnuplot "${path}/"*.gnuplot
done

exit 0