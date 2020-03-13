#!/bin/sh
######
# Copyright (c) 2020 Dr Kaifu Chen lab at the Houston Methodist Research Institute
# Author: Guangyu Wang and Shuo Zhang
# Contact:
#
######

## usage function tells users how to run stripeDiff.sh
usage() {
    echo "Usage: $0 [options]* <-a matrixA> <-b matrixB>"
    echo "Required arguments:"
    echo "  <-a matrixA> defines the input matrix A"
    echo "  <-b matrixB> defines the input matrix B"
    echo "	<-n name> defines two sample names and chromosome. e.g., wt,mutant,chr1"
    echo ""
    echo "Optional arguments:"
    echo "	[-l length] defines row and column number of split submatrix. Defult: 200"
    echo "	[-o outDir] defines the path to output files. It must be set if the two input matrices are not under the same directory"
    echo "	[-s srcDir] defines the path to stripeCalling scripts. Defualt: current working directory"
    echo ""
    exit $1
}

# get path for scripts
srcDir=$(dirname $0)

## Set default options
# default output, can be set via -o
outDir=""
# default row and column number of split submatrix is 200
length=200


## Parse command line arguments
while getopts "a:b:l:n:o:h:" opt
do
	case $opt in
		a) matrixA=$OPTARG ;;
		b) matrixB=$OPTARG ;;
		l) length=$OPTARG ;;
		n) name=$OPTARG ;;
		o) outDir=$OPTARG ;;
        h) usage 0 ;;
        [?]) usage 1 ;;
	esac
done
echo $srcDir

### Check command line arguments
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
matrixADir=$(dirname $matrixA)
matrixABase=$(basename $matrixA)
matrixBDir=$(dirname $matrixB)
matrixBBase=$(basename $matrixB)
if [ ! -z "$outDir" ]
then
	if [ ! -d $outdir ]
	then
		echo "*** error: provided output directory doesn't exist ***"
		usage 1
	fi
else
	if [ "$matrixADir" == "$matrixBDir" ]
	then
		outDir=$matrixADir
	else
		echo "*** error: the two input matrices are not under the same directory ***"
		echo "*** error: output directory must be provided ***"
		usage 1
	fi
fi
echo "The output directory is: "
echo $outDir

# check matrix name and chromsome are provided
if [ -z "$name" ]
then
	echo "*** error: matrix name and chromsome must be provided ***"
	usage 1
else
	aliasA=$(awk -F"," '{print $1}' <<< $name)
	aliasB=$(awk -F"," '{print $2}' <<< $name)
	chrom=$(awk -F"," '{print $3}' <<< $name)
fi
 ### The end of checking command line arguments


### Split matrices A and B
<<COMMENT
if [ file $matrixA | grep -q "gzip compressed data" ]
then
	gunzip $matrixA
	$matrixA =
fi
COMMENT

echo
echo "Spliting $matrixA ......"
python ${srcDir}/matrixUtils.py matrixSplit $matrixA -l $length -o ${outDir}/${aliasA}_${chrom}_splitMatrix
if [ $? != 0 ]
then
	echo "*** error: Split $matrixA failed ***"
	usage 1
else
	echo "$matrixA is split successfully!"
fi

echo "Spliting $matrixB ......"
python ${srcDir}/matrixUtils.py matrixSplit $matrixB -l $length -o ${outDir}/${aliasB}_${chrom}_splitMatrix
if [ $? != 0 ]
then
	echo "*** error: Split $matrixB failed ***"
	usage 1
else
	echo "$matrixB is split successfully!"
fi
### End of spliting matrices A and B


### Call differential stripes
cd $outDir
# The stripes directory contains called differential stripes
[ -d stripes ] && rm -r stripes
mkdir stripes
# get subchr comparisons
echo ""
echo "Extracting subchr comparisons ......"
python ${srcDir}/matrixUtils.py getSubchrComparison ${aliasA}_splitMatrix ${aliasB}_splitMatrix subchrComparison.txt
if [ $? != 0 ]
then
	echo "*** error: Extracting subchr comparisons is failed ***"
	usage 1
fi
sort -k1 -n subchrComparison.txt > sorted_subchrComparison.txt
echo "Extracting subchr comparisons is done!"
echo ""

while read -r comparison subchrA parameterA subchrB parameterB
do
	stripeOutDir="./stripes/${aliasA}_${aliasB}.${comparison}"
	echo "Calling stripes for ${stripeOutDir} ......"
	Rscript ${srcDir}/stripeDiffCalling.R -f ${subchrA},${subchrB} -p ${parameterA},${parameterB} -o ${stripeOutDir}
	if [ $? != 0 ]
	then
		echo "*** error: Calling stripes for ${stripeOutDir} is failed  ***"
		usage 1
	else
		echo "Calling stripes for ${stripeOutDir} is Done!"
	fi
	echo
done < sorted_subchrComparison.txt
### End of calling differential stripes


### combine called differential stripes
python ${srcDir}/matrixUtils.py combineStripe stripes sorted_subchrComparison.txt $name
