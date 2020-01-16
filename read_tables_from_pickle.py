#!/usr/bin/env python


from file_table_functions import *
import argparse as arg
from collections import OrderedDict
import pickle
import os, copy
import numpy as np

if __name__ == "__main__":
 
    parser = arg.ArgumentParser()
    parser.add_argument('-i', '--input_dir', type=str)
    parser.add_argument('-of', '--out_file', type=str)

    args = parser.parse_args()

    input_dir = args.input_dir
    out_file_name = args.out_file

    pickle_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if (os.path.isfile(os.path.join(input_dir, f)) and ('.pkl' in f))]

    reco_tmp_array = OrderedDict()
    reco_array = OrderedDict()

    reference = None

    for i, pick in enumerate(pickle_files):
        tmp_array = None
        with open(pick, "rb") as f:
            tmp_array = pickle.load(f)
            reference = pickle.load(f)
            for sample in tmp_array:
                if sample not in reco_tmp_array:
                    reco_tmp_array[sample] = OrderedDict()
                for tree in tmp_array[sample]:
                    if tree not in reco_tmp_array[sample]:
                        reco_tmp_array[sample][tree] = []
                    reco_tmp_array[sample][tree].append(list(tmp_array[sample][tree]))
            f.close()

    for sample in reco_tmp_array:
         reco_array[sample] = OrderedDict()

         for tree in reco_tmp_array[sample]:
             for i, num_list in enumerate(reco_tmp_array[sample][tree]):
                 if i == 0:
                     reco_array[sample][tree] = num_list
                 else:
                     for j, num in enumerate(num_list):
                         reco_array[sample][tree][j] += num

    write_table(reco_array, reference, out_file_name)

