

### combine called differential stripes
# p_cutoff: output only differential stripes whose P-value is <= p_cutoff
# indir: a directory that contains called differential stripes
# outdir: the directory that stores combined stripes
python ${StripeDiff_dir}/src/combine_stripes.py $p_cutoff $indir $outdir
