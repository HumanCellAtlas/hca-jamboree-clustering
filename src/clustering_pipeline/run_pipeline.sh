#!/bin/bash

usage(){
	echo  -e "Usage: $0 [-h] -i -o \n"
	echo "This script takes as input the location of the expression matrix"
	echo "in sparse matrix format and runs all the clustering procedures."
	echo "arguments:"
	echo " -i --inputdir     Input expression matrices location"
	echo " -o --outputdir    Location for output files"
	echo " -h --help         how help message and exit"
	exit
}

while [[ $# -gt 0 ]]
do
key="$1"
case $key in
	-i|--inputdir)
	INPUTDIR="$2"
	shift
	;;
	-o|--outputdir)
	OUTPUTDIR="$2"
	shift
	;;
	-h|--help)
	usage
	exit
	;;
esac
shift
done


if [ -z ${INPUTDIR+x} ]; then echo "Input directory not set"; usage; fi
if [ -z ${OUTPUTDIR+x} ]; then echo "Output location not set"; usage; fi


for mat in ${INPUTDIR}/*mtx; do
	echo $mat
	matname=$(basename "$mat" .mtx)
	Rscript --default-packages=stats,graphics,grDevices,utils,datasets,base,methods slm_clustering.R $mat $OUTPUTDIR $matname & Rscript --default-packages=stats,graphics,grDevices,utils,datasets,base,methods tsne_clustering.R $mat $OUTPUTDIR $matname
done
