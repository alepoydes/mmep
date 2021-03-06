#!/bin/bash

#BASE="$@"
#DESC=`echo $BASE | tr -dc '[:alnum:]\n\r'`

do_skipping=true
if [ "$1" = "-S" ]; then do_skipping=false; shift 1; fi

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

trap "{ rm -f ${lock} ; exit 255; }" EXIT

mkdir -p "${PREFIX}"
touch "${lock}"

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
	if [ "$do_skipping" = true ] && [ "$(ls -A ${path} 2>/dev/null)" ]; then
		echo -e "\e[93mSkipping\e[0m non empty $path"
	else
		echo -e "\e[92mPopulating\e[0m ${path}"
		mkdir -p "${path}"
		$@ -D "${path}" || exit 1
	fi
	gnuplot "${path}/"*.gnuplot
done < ${PARAM}

exit 0
#rm "${lock}"