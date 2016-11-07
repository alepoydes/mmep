#!/bin/bash

PREFIX="tmp/"
PARAM=$1
shift 1

DESC=`echo $@ | tr -dc '[:alnum:]\n\r'`
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
	path=$PREFIX
	for comp in "${components[@]}"; do
		IFS='=' read var val <<< "${comp}"
		echo export "$var=$val"
		export "$var=$val"
		path="${path}$var$val"
	done
	path="${path}/"
	if [ "$(ls -A ${path} 2>/dev/null)" ]; then
		echo -e "\e[93mSkipping\e[0m non empty $path"
	else
		echo -e "\e[92mPopulating\e[0m " $path
		mkdir -p "${path}"
		$@ -D "${path}" || exit 1
		gnuplot "${path}/"*.gnuplot
	fi
done < ${PARAM}

exit 0
#rm "${lock}"