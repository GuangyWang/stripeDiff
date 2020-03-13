#!/usr/bin/env python
import numpy as np
import pandas as pd
import os
import sys
import argparse


def file_split(file, resolution, outfile, length, print_subcontact, aliases):
    '''
    divide whole chromosome into different sub chromosome files and calculate cutoffs
    :param file: input chromosome file
    :param outfile: output path for sub chromosome files
    :param length: length of sub chromosome
    :param print_subcontact: if print_subcontact==1, not save sub chromosome files, if print_subcontact!=1, save sub
    chromosome files
    :return: up limit and down limit for this chromosome
    '''

    base = os.path.basename(file)
    if aliases != "None":
        filename = aliases
    else:
        filename = os.path.splitext(base)[0]
    contact_map = readMatrix(file, resolution)
    r, c = contact_map.shape
    num = int(r/length)
    mean_list = []
    files = []
    flag = 0
    if (num*length+length/2) > c:
        flag = 1
    if flag == 1:
        for i in range(num):
            sub_contact = contact_map[(i * length):((i + 1) * length), (i * length):((i + 1) * length)]
            out = os.path.join(outfile, filename + '.' + str(i * length) + '.' + str((i + 1) * length) + '.subchr')
            files.append(out)
            mean_list.append(np.mean(sub_contact))
            if print_subcontact != 1:
                np.savetxt(out, sub_contact, delimiter='\t')
            if i != (num-1):
                sub_contact2 = contact_map[int(i * length+length/2):int((i + 1) * length+length/2),
                               int(i * length+length/2):int((i + 1) * length+length/2)]
                out2 = os.path.join(outfile,
                                   filename + '.' + str(int(i * length+length/2)) + '.' + str(int((i + 1) * length + length/2)) + '.subchr')
                files.append(out2)
                mean_list.append(np.mean(sub_contact2))
                if print_subcontact != 1:
                    np.savetxt(out2, sub_contact2, delimiter='\t')
        sub_contact = contact_map[(num * length):c, (num * length):c]
        out = os.path.join(outfile, filename + '.' + str(num * length) + '.' + str(c) + '.subchr')
        files.append(out)
        mean_list.append(np.mean(sub_contact))
        if print_subcontact != 1:
            np.savetxt(out, sub_contact, delimiter='\t')
        sub_contact2 = contact_map[int((num-1) * length + length / 2):int(((num-1) + 1) * length + length / 2),
                       int((num-1) * length + length / 2):int(((num-1) + 1) * length + length / 2)]
        out2 = os.path.join(outfile, filename + '.' + str(int((num-1) * length + length / 2)) + '.' + str(int(((num-1) + 1) * length + length / 2)) + '.subchr')
        files.append(out2)
        mean_list.append(np.mean(sub_contact2))
        if print_subcontact != 1:
            np.savetxt(out2, sub_contact2, delimiter='\t')
    else:
        for i in range(num):
            sub_contact = contact_map[(i * length):((i + 1) * length), (i * length):((i + 1) * length)]
            out = os.path.join(outfile, filename + '.' + str(i * length) + '.' + str((i + 1) * length) + '.subchr')
            files.append(out)
            mean_list.append(np.mean(sub_contact))
            if print_subcontact != 1:
                np.savetxt(out, sub_contact, delimiter='\t')
            sub_contact2 = contact_map[int(i * length+length/2):int((i + 1) * length+length/2),
                           int(i * length+length/2):int((i + 1) * length+length/2)]
            out2 = os.path.join(outfile,
                               filename + '.' + str(int(i * length+length/2)) + '.' + str(int((i + 1) * length + length/2)) + '.subchr')
            files.append(out2)
            mean_list.append(np.mean(sub_contact2))
            if print_subcontact != 1:
                np.savetxt(out2, sub_contact2, delimiter='\t')
        sub_contact = contact_map[(num * length):c, (num * length):c]
        out = os.path.join(outfile, filename + '.' + str(num * length) + '.' + str(c) + '.subchr')
        files.append(out)
        mean_list.append(np.mean(sub_contact))
        if print_subcontact != 1:
            np.savetxt(out, sub_contact, delimiter='\t')
        sub_contact2 = contact_map[int(num * length + length / 2):int((num + 1) * length + length / 2),
                       int(num * length + length / 2):int((num + 1) * length + length / 2)]
        out2 = os.path.join(outfile,
                            filename + '.' + str(int(num * length + length / 2)) + '.' + str(int(
                                (num + 1) * length + length / 2)) + '.subchr')
        files.append(out2)
        mean_list.append(np.mean(sub_contact2))
        if print_subcontact != 1:
            np.savetxt(out2, sub_contact2, delimiter='\t')
    mean_list = np.array(mean_list)
    mean_list = mean_list[~np.isnan(mean_list)]
    cutoff_up = np.median(np.array(mean_list) / 0.35)
    cutoff_down = np.median(np.array(mean_list) / 2)
    files.sort()
    outputParameter = os.path.join(outfile, 'parameter.txt')
    nf = open(outputParameter, 'w')
    nf.write(str(cutoff_up)+'\t'+str(cutoff_down) + '\n')

    return cutoff_up, cutoff_down, files


def printHelp():
    print('\nmatrixSplit version 0.0.2')
    print('For help information for each function, try:\npython3 matrixSplit.py <function> -h')
    print('\nFunctions:')
    print('\tmatrixSplit:\n\t\t split Hi-C matrix\n')
    print('')


def matrixSplit(command='matrixSplit'):
    '''
    matrixSplit
    '''

    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython3 matrixSplit.py matrixSplit <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python3 matrixSplit.py matrixSplit -h\n")
        return 0

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython3 matrixSplit.py matrixSplit <contact_map_file_paths> "
                                           "[optional arguments]\n\n", description='')

    parser.add_argument('command', default=None,
                        help="set as 'matrixSplit' to split Hi-C matrix.")

    parser.add_argument('contact_map', help="Hi-C contact map.")

    parser.add_argument('-o', '--output', dest="output", default='None',
                        help="path to output files")

    parser.add_argument('-l', '--length', dest="length", default=200,
                        help="format of output files (xlsx or csv)")

    parser.add_argument('-p', '--print_subcontact', dest="print_subcontact", default=0,
                        help="format of output files (xlsx or csv)")

    parser.add_argument('--contact_maps_aliases', dest="aliases", default='None',
                        help="aliases for the contact map.")

    parser.add_argument('--resolution', dest="resolution", default='None', help="resolution for Hi-C contact map. If not provided, "
                                                                                "it will be inferred from verbose matrix")
    args = parser.parse_args()

    if args.output == 'None':
        outputPath = os.path.join(os.getcwd(), 'splitFiles')
    else:
        outputPath = args.output

    if os.path.exists(outputPath):
        print('Output path already exits')
    else:
        os.mkdir(outputPath)
        if os.path.exists(args.contact_map):
            file_split(args.contact_map, args.resolution, args.output, int(args.length), args.print_subcontact, args.aliases)
        else:
            if os.path.exists(os.path.join(os.getcwd(), args.contact_map)):
                file_split(os.path.exists(os.path.join(os.getcwd(), args.contact_map)), args.resolution,
                           outputPath, int(args.length), args.print_subcontact, args.aliases)
            else:
                print('Input file does not exit')


def readMatrix(file, resolution):
    fh = open(file)
    first_line = fh.readline()
    first_line = first_line.strip().split()
    if len(first_line) > 3: # It is a dense matrix
        fh.close()
        return np.loadtxt(file)
    else:
        max_site = 0
        inferred_resolution = 0
        for line in fh:
            (start, end, interaction) = line.strip().split()
            start = int(start)
            end = int(end)
            max_site = max(max_site, start, end)
            distance = abs(start - end)
            if inferred_resolution == 0:
                inferred_resolution = distance
            elif distance > 0 and inferred_resolution > distance: # resolution is the smallest distance between two interacting loci
                inferred_resolution = distance
        fh.close()

        # if resolution is not provided, inferred resolution will be used
        if resolution == "None":
            resolution = inferred_resolution
        # create a all-zero matrix
        dim = max_site // resolution + 1
        res = np.zeros((dim, dim))
        with open(file) as fh:
            for line in fh:
                (start, end, interaction) = line.strip().split()
                start = int(start)
                end = int(end)
                interaction = float(interaction)
                row_index = start // resolution
                col_index = end // resolution
                res[row_index][col_index] = interaction
                res[col_index][row_index] = interaction
        return res


def getSubchrComparison(command="getSubchrComparison"):
    '''
    getSubchrComparison
    '''

    if (len(sys.argv) < 4) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least three parameters need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython3 matrixUtil.py matrixSplit <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python3 matrixSplit.py matrixSplit -h\n")
        return 0

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython3 matrixSplit.py matrixSplit <contact_map_file_paths> "
                                           "[optional arguments]\n\n", description='')

    parser.add_argument('command', default='None',
                        help="set as 'getSubchrComparison' to get comparisons of two sub contact map.")

    parser.add_argument('matrixA_path', help="path to matrix A subchr")
    
    parser.add_argument('matrixB_path', help="path to matrix B subchr")
    
    parser.add_argument('outfile', help="output file")
    
    args = parser.parse_args()
    # check if the input path exists
    if not os.path.exists(args.matrixA_path):
        sys.exit("The path " + args.matrixA_path + " doesn't exist")
    if not os.path.exists(args.matrixB_path):
        sys.exit("The path " + args.matrixB_path + " doesn't exist")

    # get parameters for matrix A and B
    parameterFileA = args.matrixA_path + '/parameter.txt'
    parameterA = getParameter(parameterFileA)
    parameterFileB = args.matrixB_path + '/parameter.txt'
    parameterB = getParameter(parameterFileB)

    # get subchr
    outfh = open(args.outfile, 'w')
    subchr_A = [f for f in os.listdir(args.matrixA_path) if f.endswith('subchr')]
    subchr_B = [f for f in os.listdir(args.matrixB_path) if f.endswith('subchr')]
    subchr_A.sort()
    subchr_B.sort()
    if len(subchr_A) != len(subchr_B):
        print("The two input matrices do not have the same size")
        exit(1)
    for i in range(len(subchr_A)):
        file_A = subchr_A[i]
        file_B = subchr_B[i]
        temp = file_A.split('.')
        positionA = temp[-3] + '.' + temp[-2]
        temp = file_B.split('.')
        positionB = temp[-3] + '.' + temp[-2]
        if positionA != positionB:
            print(positionA)
            print("The two input matrices do not have the same size")
            exit(1)
        outfh.write(positionA + '\t' + args.matrixA_path + '/' + file_A + '\t' + parameterA + '\t' + \
            args.matrixB_path + '/' + file_B + '\t' + parameterB + '\n')
    outfh.close()


def combineStripe(command="combineStripe"):
    '''
    combineStripe
    '''

    if (len(sys.argv) < 4) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least three parameters need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython3 matrixUtil.py matrixSplit <contact_map_file_paths> [optional arguments]\n\n"
              "for more help, please try: python3 matrixSplit.py matrixSplit -h\n")
        return 0

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython3 matrixSplit.py matrixSplit <contact_map_file_paths> "
                                           "[optional arguments]\n\n", description='')

    parser.add_argument('command', default='None',
                        help="combine differential stripes")

    parser.add_argument('data_path', help="path to differential stripes")
    parser.add_argument('comparison', help="get comparison position")
    
    parser.add_argument('name', help="two sample names and chromosome. e.g., wt,mutant,chr1")
    
    args = parser.parse_args()
    (sample1, sample2, chrom) = args.name.split(',')
    # check if the input path exists
    if not os.path.exists(args.data_path):
        sys.exit("The path " + args.data_path + " doesn't exist")

    if not args.data_path.endswith('/'):
        args.data_path += '/'
    
    # write output header
    in_sample1 = "in_" + sample1 + "_not_" + sample2 + "_stripes.txt"
    outfh = open(in_sample1, 'w')
    outfh.write("chrom\t" + "upPeak.loc" + '\t' + "downPeak.loc" + '\t' + "leftEdge" + '\t'+ "rightEdge" + '\t' + "upPeak.sample1" + \
        '\t' + "downPeak.sample1" + '\t' + "logFoldChange.sample1" + '\t' + "strap.pValue.sample1" + '\t' + "upPeak.sample2" + \
        '\t' + "downPeak.sample2" + '\t' + "logFoldChange.sample2" + '\t' + "strap.pValue.sample2" + '\t' + "diffStrap.pValue" + \
        '\t' + "direction" + '\n')
    outfh.close()

    in_sample2 = "in_" + sample2 + "_not_" + sample1 + "_stripes.txt"
    outfh = open(in_sample2, 'w')
    outfh.write("chrom\t" + "upPeak.loc" + '\t' + "downPeak.loc" + '\t' + "leftEdge" + '\t'+ "rightEdge" + '\t' + "upPeak.sample1" + \
        '\t' + "downPeak.sample1" + '\t' + "logFoldChange.sample1" + '\t' + "strap.pValue.sample1" + '\t' + "upPeak.sample2" + \
        '\t' + "downPeak.sample2" + '\t' + "logFoldChange.sample2" + '\t' + "strap.pValue.sample2" + '\t' + "diffStrap.pValue" + \
        '\t' + "direction" + '\n')
    outfh.close()

    # combine stripe
    mylist = readComparison(args.comparison)
    for position in mylist:
        start = int(position.split('.')[0])
        infile1 = args.data_path + sample1 + '_' + sample2 + '.' + position + "/1.txt"
        infile2 = args.data_path + sample1 + '_' + sample2 + '.' + position + "/2.txt"
        if (not os.path.exists(infile1)) or (not os.path.exists(infile2)):
            continue
        reformat(infile1, start, chrom, in_sample1)
        reformat(infile2, start, chrom, in_sample2)



def readComparison(infile):
    mylist = []
    with open(infile) as f:
        for line in f:
            position = line.strip().split()
            mylist.append(position[0])
    return mylist


def reformat(infile, start, chrom, outfile):
    outfh = open(outfile, 'a+')
    with open(infile) as f:
        for row in f:
            row = row.strip().split()
            if "upPeak.loc" in row[0] or "NA" in row:
                continue
            line = chrom + '\t'
            for i in range(4):
                line += str(int(row[i]) + start) + '\t'
                for i in range(4, 14):
                    line += str(row[i])
            outfh.write(line + '\n')
    outfh.close()


def getParameter(infile):
	with open(infile) as f:
		line = f.readline()
		line = line.strip().split()
		return line[0]


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == 'matrixSplit':
            matrixSplit(command='matrixSplit')
        elif sys.argv[1] == "getSubchrComparison":
            getSubchrComparison(command="getSubchrComparison")
        elif sys.argv[1] == "combineStripe":
            combineStripe(command="combineStripe")
        else:
            printHelp()
    else:
        print('\nmatrixSplit version 0.0.2')
        print('For a list of functions in matrixSplit, please try:\npython3 matrixSplit.py -h')
        print('')

