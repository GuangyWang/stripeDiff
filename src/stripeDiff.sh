#!/bin/sh
######
# Copyright (c) 2020 Dr Kaifu Chen lab at the Houston Methodist Research Institute
# Author: Guangyu Wang and Shuo Zhang
# Contact: gwang2@houstonmethodist.org and szhang3@houstonmethodist.org
#
######

### usage function tells users how to run stripeDiff.sh
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

# Get absolute path for scripts and check if required scripts exist
srcDir=$(dirname $0)
cd $srcDir
srcDir=$PWD
cd -
# check if matrixUtils.py exists
if [ ! -e $srcDir/matrixUtils.py ]
then
	echo "matrixUtils.py is not exist in $srcDir, please put all scripts in the directory containing stripeDiff.sh"
	usage 1
fi
# check if the strieDiff.R exists
if [ ! -e $srcDir/stripeDiffCalling.R ]
then
	echo "stripeDiffCalling.R is not exist in $srcDir, please put all scripts in the directory containing stripeDiff.sh\n"
	usage 1
fi


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
	if [ ! -d $outDir ]
	then
		mkdir $outDir
		if [ $? != 0 ]
		then
			echo "*** error: provided output directory doesn't exist and cannot be created ***"
			usage 1
		fi
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

# Check if matrix name and chromsome are provided
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

echo "Spliting $matrixA ......"
python3 ${srcDir}/matrixUtils.py matrixSplit $matrixA -l $length -o ${outDir}/${aliasA}_${chrom}_splitMatrix
if [ $? != 0 ]
then
	echo "*** error: Split $matrixA failed ***"
	usage 1
else
	echo "$matrixA is split successfully!"
fi

echo "Spliting $matrixB ......"
python3 ${srcDir}/matrixUtils.py matrixSplit $matrixB -l $length -o ${outDir}/${aliasB}_${chrom}_splitMatrix
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
[ -d ${chrom}_stripes ] && rm -r ${chrom}_stripes
mkdir ${chrom}_stripes
# get subchr comparisons
echo ""
echo "Extracting subchr comparisons ......"
python3 ${srcDir}/matrixUtils.py getSubchrComparison ${aliasA}_${chrom}_splitMatrix ${aliasB}_${chrom}_splitMatrix subchrComparison.txt
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
	stripeOutDir="./${chrom}_stripes/${aliasA}_${aliasB}.${comparison}"
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
# output: in_${aliasA}_not_${aliasB}_${chrom}_stripes.txt and in_${aliasA}_not_${aliasB}_${chrom}_stripes.txt
# get resolution
resolution=10000
python3 ${srcDir}/matrixUtils.py combineStripe ${chrom}_stripes sorted_subchrComparison.txt $resolution $name
# sort differential stripes based on estimated position
sort -k4,5 -n in_${aliasA}_not_${aliasB}_${chrom}_stripes.txt > sorted_in_${aliasA}_not_${aliasB}_${chrom}_stripes.txt
sort -k4,5 -n in_${aliasB}_not_${aliasA}_${chrom}_stripes.txt > sorted_in_${aliasB}_not_${aliasA}_${chrom}_stripes.txt

# remove duplications
infile=sorted_in_${aliasA}_not_${aliasB}_${chrom}_stripes.txt
outfile=deduplicated_in_${aliasA}_not_${aliasB}_${chrom}_stripes.txt
python3 ${srcDir}/matrixUtils.py deduplicate $infile $outfile
infile=sorted_in_${aliasB}_not_${aliasA}_${chrom}_stripes.txt
outfile=deduplicated_in_${aliasB}_not_${aliasA}_${chrom}_stripes.txt
python3 ${srcDir}/matrixUtils.py deduplicate $infile $outfile
### End of combining called differential stripes
