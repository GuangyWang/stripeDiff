#!/bin/sh
######
# Copyright (c) 2020 Dr. Kaifu Chen lab
# Author: Guangyu Wang and Shuo Zhang
# Contact: gwang2@houstonmethodist.org
#
######


### usage function tells users how to run stripeDiff.sh
usage() {
	echo "***************** How to use stripeDiff *****************"
    echo "Usage: sh /path/to/stripeDiff.sh [options]* <-a matrixA> <-b matrixB>"
    echo "Required arguments:"
    echo " 	<-a matrixA> defines the input matrix A"
    echo " 	<-b matrixB> defines the input matrix B"
    echo "	<-n name> defines two sample names and chromosome. e.g., wt,mutant,chr1"
    echo ""
    echo "Optional arguments:"
    echo "	[-l length] defines row and column number of split submatrix. Defult: 300"
    echo "	[-o outDir] defines the path to output files. It must be set if the two input matrices are not under the same directory"
    echo "	[-r resolution] defines the bin size in base pair. If it is not provided, the bin number of identified stripes will be output"
    echo "	[-f] runs in a faster mode without estimating stripe length\n\n"
    echo "For detailed informatiom, please feel free to refer: https://github.com/GuangyWang/stripeDiff"
    echo "*********************************************************"
    echo ""
    exit $1
}

### Set default options
# default output, can be set via -o
outDir=""
# default row and column number of split submatrix is 200. It can be set via -l
length=300
# default bin size, can be set via -r
resolution=""
# estimate stripe length. To turn off, use -f
estimateLen=1
### The end of setting default options


### Parse command line arguments
while getopts "a:b:l:n:o:fhr:" opt
do
	case $opt in
		a) matrixA=$OPTARG ;;
		b) matrixB=$OPTARG ;;
		l) length=$OPTARG ;;
		n) name=$OPTARG ;;
		o) outDir=$OPTARG ;;
		r) resolution=$OPTARG ;;
		f) estimateLen=0 ;;
        h) usage 0 ;;
        [?]) usage 1 ;;
	esac
done


date
echo "The following command is running:"
echo "	sh $0 $*\n"

# Get absolute path for scripts and check if required scripts exist
temp=$PWD
srcDir=$(dirname $0)
cd $srcDir
srcDir=$PWD
cd $temp


### check if all required scripts exist in desired directory
# check if matrixUtils.py exists
if [ ! -e "$srcDir"/matrixUtils.py ]
then
	echo "*** error: matrixUtils.py is not exist in $srcDir, please put all scripts in the directory containing stripeDiff.sh ***\n"
	usage 1
fi
# check if the stripe_detect.R exists
if [ ! -e "$srcDir"/stripe_detect.R ]
then
	echo "*** error: stripe_detect.R is not exist in $srcDir, please put all scripts in the directory containing stripeDiff.sh ***\n"
	usage 1
fi
# check if the differential_score.R exists
if [ ! -e "$srcDir"/differential_score.R ]
then
	echo "*** error: differential_score.R is not exist in $srcDir, please put all scripts in the directory containing stripeDiff.sh ***\n"
	usage 1
fi
### End of checking requried scripts


### Check command line arguments
# Check if matrix A is provided
if [ -z "$matrixA" ]
then
	echo "*** error: input matrix A must be provide ***"
	usage 1
else
	if [ ! -e "$matrixA" ]
	then
		echo "*** error: provided matrix A doesn't exist ***"
		usage 1
	fi
fi

# Check if matrix B is provided
if [ -z "$matrixB" ]
then
	echo "*** error: input matrix B must be provide ***"
	usage 1
else
	if [ ! -e "$matrixA" ]
	then
		echo "*** error: provided matrix B doesn't exist ***"
		usage 1
	fi
fi

# Check if matrix A and matrix B are stored in a same directory. If not, output directory must be provided
matrixADir=$(dirname $matrixA)
matrixABase=$(basename $matrixA)
matrixBDir=$(dirname $matrixB)
matrixBBase=$(basename $matrixB)
if [ ! -z "$outDir" ]
then
	if [ ! -d "$outDir" ]
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
	echo "*** error: matrix name and chromosome must be provided ***"
	usage 1
else
	aliasA=$(awk -F"," '{print $1}' <<< $name)
	aliasB=$(awk -F"," '{print $2}' <<< $name)
	chrom=$(awk -F"," '{print $3}' <<< $name)
	if [ -z "$aliasA" -o -z "$aliasB" -o -z "$chrom" ]
	then
		echo "*** error: matrix name and chromosome should be separated by ',' and without any space between them ***"
		echo "***        e.g., mt,mutant,chr1 ***"
		usage 1
	fi
fi

# check if the bin size is provided. If not, set bin size to 1
if [ -z "$resolution" ]
then
	resolution=1
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
# determine if estimating stripe length is disabled
rCode=""
if [ "$estimateLen" == 1 ]
then
	rCode="stripe_detect.R"
else
	rCode="stripe_detect_no_length.R"
fi
rCode="test_noLength.R"
# The stripes directory contains called differential stripes
[ -d "${chrom}_stripes" ] && rm -r ${chrom}_stripes
mkdir ${chrom}_stripes
# get subchr comparisons
echo ""
echo "Extracting subchr comparisons ......"
python ${srcDir}/matrixUtils.py getSubchrComparison ${aliasA}_${chrom}_splitMatrix ${aliasB}_${chrom}_splitMatrix subchrComparison_${chrom}.txt
if [ $? != 0 ]
then
	echo "*** error: Extracting subchr comparisons is failed ***"
	usage 1
fi
sort -k1 -n subchrComparison_${chrom}.txt > sorted_subchrComparison_${chrom}.txt
echo "Extracting subchr comparisons is done!"
echo ""

while read -r comparison subchrA parameterA subchrB parameterB
do
	stripeOutDir="./${chrom}_stripes/${aliasA}_${aliasB}.${comparison}"
	echo "Calling stripes for ${stripeOutDir} ......"
	echo "Running command: Rscript ${srcDir}/${rCode} -f ${subchrA},${subchrB} -p ${parameterA},${parameterB} -o ${stripeOutDir}"
	Rscript ${srcDir}/${rCode} -f ${subchrA},${subchrB} -p ${parameterA},${parameterB} -o ${stripeOutDir}
	if [ $? != 0 ]
	then
		echo "*** error: Calling stripes for ${stripeOutDir} is failed  ***"
	else
		echo "Calling stripes for ${stripeOutDir} is Done!"
	fi
	echo
done < sorted_subchrComparison_${chrom}.txt
### End of calling differential stripes


### combine called differential stripes
# output: in_${aliasA}_not_${aliasB}_${chrom}_stripes.txt and in_${aliasA}_not_${aliasB}_${chrom}_stripes.txt
python ${srcDir}/matrixUtils.py combineStripe ${chrom}_stripes sorted_subchrComparison_${chrom}.txt $resolution $name --estimateLen $estimateLen
# sort differential stripes based on estimated position
sort -k4,5 -n in_${aliasA}_not_${aliasB}_${chrom}_stripes.txt > sorted_in_${aliasA}_not_${aliasB}_${chrom}_stripes.txt
sort -k4,5 -n in_${aliasB}_not_${aliasA}_${chrom}_stripes.txt > sorted_in_${aliasB}_not_${aliasA}_${chrom}_stripes.txt
### End of combining


### calculate P values for each differential stripe
infile="sorted_in_${aliasA}_not_${aliasB}_${chrom}_stripes.txt"
outfile="in_${aliasA}_not_${aliasB}_${chrom}_stripes_with_p.txt"
Rscript ${srcDir}/differential_score.R -f $infile -o $outfile

infile="sorted_in_${aliasB}_not_${aliasA}_${chrom}_stripes.txt"
outfile="in_${aliasB}_not_${aliasA}_${chrom}_stripes_with_p.txt"
Rscript ${srcDir}/differential_score.R -f $infile -o $outfile
### End of calculating P values


# remove duplications
infile="in_${aliasA}_not_${aliasB}_${chrom}_stripes_with_p.txt"
outfile="deduplicated_in_${aliasA}_not_${aliasB}_${chrom}_stripes.txt"
python ${srcDir}/matrixUtils.py deduplicate $infile $outfile --estimateLen $estimateLen

infile="in_${aliasB}_not_${aliasA}_${chrom}_stripes_with_p.txt"
outfile="deduplicated_in_${aliasB}_not_${aliasA}_${chrom}_stripes.txt"
python ${srcDir}/matrixUtils.py deduplicate $infile $outfile --estimateLen $estimateLen
### End of combining called differential stripes


echo "Calling differential stripes is done!"
date
