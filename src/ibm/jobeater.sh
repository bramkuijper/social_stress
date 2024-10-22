#!/usr/bin/env bash

# jobeater.sh - execute a batch file $1 over $2 cores.

if [ "$#" -ne 2 ];
	     then echo "1. provide file name and 2. number of cores you want to
		     be using."
		          exit
fi

cat $1 | tr '\n' '\0' | xargs -0 -P$2 -n1 sh -c
