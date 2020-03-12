#!/bin/sh
######
# Copyright (c) 2020 Dr Kaifu Chen lab at the Houston Methodist Research Institute



######

## usage function tells users how to run stripeDiff.sh
usage() {
    echo "Usage: $0 [options]* <-a matrixA> <-b matrixB>"
    echo "Required arguments:"
    echo "  <-a matrixA> defines the input matrix A"
    echo "  <-b matrixB> defines the input matrix B"
    echo ""
    echo "Optional arguments:"
    echo "	[-l length] defines row and column number of split submatrix"
    echo "	[-o outDir] defines the path to output files. It must be set if the two input matrices are not under the same directory"
    exit $1
}

## Set default options
# default output, can be set via -o
outDir=""


## Parse command line arguments
while getopts "i:o:g:p:h" opt
do
	case $opt in
		a) matrixA=$OPTARG ;;
		b) matrixB=$OPTARG ;;
		l) length=$OPTARG ;;
		o) outDir=$OPTARG ;;
		g) refGenome=$OPTARG ;;
		p) refP=$OPTARG ;;
        h) usage 0;;
        [?]) usage 1;;
	esac
done


## Check command line arguments
# Check if matrix A is provided
if [ ! -e $matrixA ]
then
	echo "*** error: input matrix A must be provide or provided matrix A doesn't exist ***"
	usage 1
fi

# Check if matrix B is provided
if [ ! -e $matrixB ]
then
	echo "*** error: input matrix B must be provide or provided matrix B doesn't exist ***"
	usage 1
fi

# Check if matrix A and matrix B are stored in a same directory. If not, output directory must be provided
if [ "$outDir" != ""]
then
	if [ ! -d $outdir ]
	then
		echo "*** error: provided output directory doesn't exist ***"
		usage 1
	fi
else
	matrixADir=$(dirname matrixA)
	matrixBDir=$(dirname matrixB)
	if [ "$matrixADir" == "$matrixBDir" ]
	then
		outDir=matrixADir
	else
		echo "*** error: the two input matrices are not under the same directory ***"
		echo "*** error: output directory must be provided ***"
		usage 1
	fi
fi

## The end of checking command line arguments

## Split matrices A and B
matrixSplit_0.0.3.py matrixSplit -c $matrixA -o ${outDir}/splitMatrixA
if [ $? != 0 ]
then
	echo "*** error: Split matrix A failed ***"
	usage 1
else
	echo "Matrix A is split successfully"
fi
matrixSplit_0.0.3.py matrixSplit -c $matrixB -o ${outDir}/splitMatrixB
if [ $? != 0 ]
then
	echo "*** error: Split matrix B failed ***"
	usage 1
else
	echo "Matrix B is split successfully"
fi

## Call differential stripes


### combine called differential stripes
# p_cutoff: output only differential stripes whose P-value is <= p_cutoff
# indir: a directory that contains called differential stripes
# outdir: the directory that stores combined stripes
python ${StripeDiff_dir}/src/combine_stripes.py $p_cutoff $indir $outdir
