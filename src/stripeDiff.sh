Usage: 

Options:
	-f CHARACTER, --inputFile=CHARACTER
		Two input Hi-C contact maps (separated by comma)

	-p CHARACTER, --parameter=CHARACTER
		Two parameters for Hi-C contact maps (separated by comma)

	-o CHARACTER, --output=CHARACTER
		Output directory name

	-h, --help
		Show this help message and exit



### combine called differential stripes
# p_cutoff: output only differential stripes whose P-value is <= p_cutoff
# indir: a directory that contains called differential stripes
# outdir: the directory that stores combined stripes
python ${StripeDiff_dir}/src/combine_stripes.py $p_cutoff $indir $outdir
