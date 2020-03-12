#!/usr/bin/env python
import numpy as np
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
    print('\nmatrixSplit version 0.0.3')
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

    parser.add_argument('
    
    
    
    
    contact_map', help="Hi-C contact map.")

    parser.add_argument('-o', '--output', dest="output", default='None',
                        help="path to output files")

    parser.add_argument('-l', '--length', dest="length", default=200,
                        help="row and column number of subcontact")

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
        print('output path already exits')
    else:
        os.mkdir(outputPath)
        if os.path.exists(args.contact_map):
            file_split(args.contact_map, args.resolution, args.output, int(args.length), args.print_subcontact, args.aliases)
        else:
            if os.path.exists(os.path.join(os.getcwd(), args.contact_map)):
                file_split(os.path.exists(os.path.join(os.getcwd(), args.contact_map)), args.resolution,
                           outputPath, int(args.length), args.print_subcontact, args.aliases)
            else:
                print('input file does not exit')


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


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == 'matrixSplit':
            matrixSplit(command='matrixSplit')
        else:
            printHelp()
    else:
        print('\nmatrixSplit version 0.0.3')
        print('For a list of functions in matrixSplit, please try:\npython3 matrixSplit.py -h')
        print('')
