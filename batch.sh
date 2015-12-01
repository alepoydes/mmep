#!/bin/bash
OPTIMIZE=./mepd
DESC=$1
SUFF=$2
shift 2
while IFS='' read -r line || [[ -n "$line" ]]; do
	echo -e "\n$line" | cat ${DESC} - | ${OPTIMIZE} $@
done < ${SUFF}