#!/usr/bin/python
import sys
import os
import pandas as pd
import numpy as np

def main(argv):
	(p_cutoff, data_path, out_path) = argv
	if not os.path.exists(out_path):
		os.mkdir(out_path, mode=0o755)

	chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']
	samples = [['CH12_wt', 'CH12_CTCFZF9to11_mutant']]
	BP = "10000" #resolution

	for s in samples:
		(stage1, stage2) = s
		in_stage1 = out_path + "in_" + stage1 + "_not_" + stage2 + '_' + p_cutoff + "_stripes.txt"
		in_stage2 = out_path + "in_" + stage2 + "_not_" + stage1 + '_' + p_cutoff + "_stripes.txt"
		for c in chroms:
			dirs = getSubchrDir(data_path, stage1, stage2, c)
			print(dirs)
			for dir in dirs:
				start = int(dir.split('.')[1])
				infile1 = data_path + dir + "/1.txt"
				infile2 = data_path + dir + "/2.txt"
				if os.path.exists(infile1):
					reformat(infile1, start, c, float(p_cutoff), int(BP), in_stage1)
				if os.path.exists(infile2):
					reformat(infile2, start, c, float(p_cutoff), int(BP), in_stage2)


def getSubchrDir(indir, stage1, stage2, chrom):
	prefix = stage1 + '_' + stage2 + '_' + chrom + '_10000'
	dirs = [ dir for dir in os.listdir(indir) if dir.startswith(prefix) ]
	return dirs


def reformat(infile, start, chrom, p_cutoff, BP, outfile):
	df = pd.read_csv(infile, delimiter='\t')
	outfh = open(outfile, 'a+')
	for row in df.itertuples():
		if pd.isnull(row[1]) or row[13] > p_cutoff:
			continue
		else:
			line = ''
			for i in range(1,5):
				line += str((int(row[i]) + start)*BP) + '\t'
			for i in range(5, 15):
				line += str(row[i]) + '\t'
			line += chrom
			outfh.write(line + '\n')
	outfh.close()
	return


def parse_arguments():
    parser = argparse.ArgumentParser()
	parser.add_argument("indir", help="The directory that contains called differential straps")
    parser.add_argument("-d", "--data_path", help="The directory that contains output from strapCalling")
    parser.add_argument("--BP", type=int, default=10000, help="BP is the resolution")
    parser.add_argument("--source", type=str, default="control", choices=["control", "TADsplimer", "random"])
    parser.add_argument("--sub", default=False, action="store_true")
    parser.add_argument("--outfile", help="output file")
    return parser


if __name__ == "__main__":
	main(sys.argv[1:])

