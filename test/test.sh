#!/bin/sh

gunzip wt.txt.gz
gunzip mutant.txt.gz
sh ../src/stripeDiff.sh  -a wt.txt -b mutant.txt -n wt,mutant,chr19
gzip wt.txt
gzip mutant.txt
