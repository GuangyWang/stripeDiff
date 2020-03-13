#!/bin/sh
gunzip CH12_CTCFZF9to11_mutant_chr19_chr19_10000_verbose.txt.gz
gunzip CH12_wt_chr19_chr19_10000_verbose.txt.gz
sh /path_to_stripeDiff/stripeDiff.sh -a CH12_wt_chr19_chr19_10000_verbose.txt -b CH12_CTCFZF9to11_mutant_chr19_chr19_10000_verbose.txt -n wt,mutant,chr19
