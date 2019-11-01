#!/usr/bin/env python

"""
Creator: Erich Schmitz
Date: Feb 22, 2019
"""

import ROOT as rt
import numpy as np
import root_numpy as rnp
import numpy.lib.recfunctions as rfc
import os
from collections import OrderedDict

rt.gROOT.SetBatch()
rt.TH1.AddDirectory(rt.kFALSE)


def get_tree_info_singular(sample_, file_name_, tree_name_, variable_list_, cuts_to_apply_=None):
    """
    Same as get_tree_info_plural, but runs over a single file
    returns structured array containing the list of variables
    """
    tmp_f = rt.TFile(file_name_, 'r')
    tmp_t = tmp_f.Get(tree_name_)
    tmp_array = None
    if bool(tmp_t) and tmp_t.InheritsFrom(rt.TTree.Class()): 
        tmp_array = rnp.tree2array(tmp_t, branches=variable_list_, selection=cuts_to_apply_)
    else:
        print('tree: ' + tree_name_ + ' is not a tree, skipping')
    
    return tmp_array


def get_tree_info_singular_deprecated(sample_, file_name_, tree_names_, variable_list_, cuts_to_apply_=None):
    """
    Same as get_tree_info_plural, but runs over a single file
    returns structured array containing the list of variables
    """
    tmp_array = {}
    tmp_array[sample_] = OrderedDict()
    for tree in tree_names_:
        tmp_f = rt.TFile(file_name_, 'r')
        tmp_t = tmp_f.Get(tree)
        if bool(tmp_t) and tmp_t.InheritsFrom(rt.TTree.Class()):
            if variable_list_: 
                tmp_array[sample_][tree] = rnp.tree2array(tmp_t, branches=variable_list_, selection=cuts_to_apply_)
            else:
                tmp_array[sample_][tree] = rnp.tree2array(tmp_t, selection=cuts_to_apply_)
        else:
            print('tree: ' + tree + ' is not a tree, skipping')
    
    return tmp_array

def reduce_and_condense(file_list_of_file_lists, variable_list):
    print('for posterity')

    tree_chain = rt.TChain('deepntuplizer/tree')
    for file_list in file_list_of_file_lists:
        files = open(file_list).readlines()
        files = [f.replace('\n','') for f in files]
        for file_name in files:
            tree_chain.Add(file_name)
        branches = [b.GetName() for b in tree_chain.GetListOfBranches()]
    for branch in branches:
        if branch not in variable_list:
            tree_chain.SetBranchStatus(branch, 0)
    file_out = rt.TFile('output_condensed.root', 'recreate')
    reduced_tree = tree_chain.CloneTree()
    reduced_tree.Write()
    file_out.Close()


def reduce_singular(in_file_name_, out_file_name_, variable_list_):
    tmp_f = rt.TFile(int_file_name_, 'r')
    tmp_t = tmp_f.Get('deepntuplizer/tree')

    branches = [b.GetName() for b in tmp_t.GetListOfBranches()]        
    for branch in branches:
        if branch not in variable_list:
            tmp_t.SetBranchStatus(branch, 0)

    out_file = rt.TFile(out_file_name_, 'recreate')
    reduced_tree = tmp_t.CloneTree()
    reduced_tree.Write()
    file_out.Close()

def process_the_samples(input_sample_list_, truncate_file_ = None, tree_in_dir_ = None):
    list_of_files = OrderedDict()

    for sample, list_of_folders in input_sample_list_.items():
        print(sample)
        file_list = []
        for folder in list_of_folders:
            print('->', folder)
            if '.root' in folder:
                file_list.append(folder)
            else:
                file_list_tmp = [os.path.join(folder, f) for f in os.listdir(folder) if (os.path.isfile(os.path.join(folder, f)) and ('.root' in f))]
                file_list.append(file_list_tmp)
        need_to_concat = [True if type(x) == list else False for x in file_list]
        if np.any(need_to_concat):
            file_list = np.concatenate(file_list)
        # Get file structure, in case there is a grid of mass points
        print file_list
        f_struct_tmp = rt.TFile(file_list[0], 'r')
        tree_list = []
        if tree_in_dir_ is not None and 'SMS' not in sample:
            for tree in tree_in_dir_:
                tree_list.append(tree)
        else:
            tree_list = [tree.GetName() for tree in f_struct_tmp.GetListOfKeys() if 'EventCount' not in tree.GetName()]

        if truncate_file_ is not None:
            file_list = file_list[:truncate_file_]
        if 'SMS' in sample:
            trees_to_keep = ['SMS_500_420', 'SMS_500_490']
            #trees_to_keep = []
            tree_name_mass = [(int(mass.split('_')[1]), int(mass.split('_')[2])) for mass in tree_list]
            tree_name_mass.sort(key=lambda x: int(x[0]))
            if trees_to_keep:
                tree_list = trees_to_keep
            else:
                tree_list = ['SMS_'+str(mom) + '_' + str(child) for mom, child in tree_name_mass]

        f_struct_tmp.Close()
        list_of_files[sample] = OrderedDict([('files', file_list), ('trees', tree_list)])
    return list_of_files



def write_hists_to_file(hists_, out_file_name_):

    out_file = rt.TFile.Open(out_file_name_, "recreate")
    print('writing histograms to: ', out_file.GetName())
    out_file.cd()
    for sample in hists_:
        sample_dir = out_file.mkdir(sample)
        for tree in hists_[sample]:
            sample_dir.cd()
            tmp_dir = sample_dir.mkdir(tree)
            tmp_dir.cd()
            for hist in hists_[sample][tree].values():
                hist.Write()
    out_file.Close()
    print('finished writing')

def write_arrays_to_file(arrays_, out_file_name_):
    """
    save ndarray to .npy file
    """
    if '.npy' in out_file_name_:
        np.save(out_file_name_, arrays_)
    else:
        raise ValueError(out_file_name_.split('.')[-1] + ' is the wrong file type, please use npy')


def evaluateZbi(Nsig, Nbkg,sys):
    Nobs = rt.Double(Nsig+Nbkg)
    tau = rt.Double(1./Nbkg/(sys*sys/10000.))
    aux = rt.Double(Nbkg*tau)
    Pvalue = rt.TMath.BetaIncomplete(1./(1.+tau),Nobs,aux+1.)
    return rt.TMath.Sqrt(2.)*rt.TMath.ErfcInverse(Pvalue*2)


def write_table(table_array_, reference_, table_name_):

    out_lines = []
    out_w_lines = []
    out_lines_zbi = []
    background = OrderedDict()
    new_signal_array = OrderedDict()
    is_background = False

    for factor, sample in enumerate(table_array_):
        if factor == 0:
            out_lines.append('{:^30}'.format(' '))
        for itree, tree in enumerate(table_array_[sample]):
            row_string = '{:<15} {:<15}'.format(sample, tree)
            out_lines.append(row_string)
            row_index = out_lines.index(row_string)
            for icol, name in enumerate(reference_):
                if '_w' in name: continue
                if factor == 0 and itree == 0:
                    out_lines[factor] += '{:^30}'.format(name)
                    out_lines[row_index] += '{:^30.4f}'.format(table_array_[sample][tree][icol])
                else:
                    out_lines[row_index] += '{:^30.4f}'.format(table_array_[sample][tree][icol])
    
    for factor, sample in enumerate(table_array_):
        if 'SMS' not in sample:
            is_background = True
        if factor == 0:
            out_w_lines.append('{:^30}'.format(' '))
        for itree, tree in enumerate(table_array_[sample]):
            row_string = '{:<15} {:<15}'.format(sample, tree)
            out_w_lines.append(row_string)
            row_index = out_w_lines.index(row_string)
            for icol, name in enumerate(reference_):
                if '_w' not in name: continue
                if is_background:
                    if name not in background:
                        background[name] = 0.
                    background[name] += table_array_[sample][tree][icol] 
                if factor == 0 and itree == 0:
                    out_w_lines[factor] += '{:^30}'.format(name)
                    out_w_lines[row_index] += '{:^30.4f}'.format(table_array_[sample][tree][icol])
                else:
                    out_w_lines[row_index] += '{:^30.4f}'.format(table_array_[sample][tree][icol])
            if not is_background:
                new_signal_array[sample+'_'+tree] = table_array_[sample][tree]
        is_background = False
    print background
    print new_signal_array
    if background:
        for factor, sample in enumerate(new_signal_array):
            if factor == 0:
                out_lines_zbi.append('{:^30}'.format(' '))
            row_string = '{:^30}'.format(sample)
            out_lines_zbi.append(row_string)
            row_index = out_lines_zbi.index(row_string)
            for icol, name in enumerate(reference_):
                if '_w' not in name: continue
                val_zbi = evaluateZbi(new_signal_array[sample][icol], background[name], 10)
                if factor == 0:
                    out_lines_zbi[factor] += '{:^30}'.format(name)
                    out_lines_zbi[row_index] += '{:^30.4f}'.format(val_zbi)
                else:
                    out_lines_zbi[row_index] += '{:^30.4f}'.format(val_zbi)

    for iline in xrange(len(out_lines)):
        out_lines[iline] += '\n'
    for iline in xrange(len(out_w_lines)):
        out_w_lines[iline] += '\n'
    for iline in xrange(len(out_lines_zbi)):
        out_lines_zbi[iline] += '\n'

    with open(table_name_, 'w') as t:
        t.writelines(out_lines)
        t.write('\n')
        t.writelines(out_w_lines)
        t.write('\n')
        t.writelines(out_lines_zbi)
        t.close()


def write_table_deperecated(table_array_, table_w_array_, table_name_):

    out_lines = []
    out_w_lines = []
    for factor, sample in enumerate(table_array_):
        if factor == 0:
            out_lines.append('{:^30}'.format(' '))
        for itree, tree in enumerate(table_array_[sample]):
            row_string = '{:<15} {:<15}'.format(sample, tree)
            out_lines.append(row_string)
            row_index = out_lines.index(row_string)
            for icol, name_value in enumerate(table_array_[sample][tree].items()):
                if factor == 0 and itree == 0:
                    out_lines[factor] += '{:^30}'.format(name_value[0])
                    out_lines[row_index] += '{:^30.4f}'.format(name_value[1])
                else:
                    out_lines[row_index] += '{:^30.4f}'.format(name_value[1])
    
    for factor, sample in enumerate(table_w_array_):
        if factor == 0:
            out_w_lines.append('{:^30}'.format(' '))
        for itree, tree in enumerate(table_w_array_[sample]):
            row_string = '{:<15} {:<15}'.format(sample, tree)
            out_w_lines.append(row_string)
            row_index = out_w_lines.index(row_string)
            for icol, name_value in enumerate(table_w_array_[sample][tree].items()):
                if factor == 0 and itree == 0:
                    out_w_lines[factor] += '{:^30}'.format(name_value[0])
                    out_w_lines[row_index] += '{:^30.4f}'.format(name_value[1])
                else:
                    out_w_lines[row_index] += '{:^30.4f}'.format(name_value[1])

    for iline in xrange(len(out_lines)):
        out_lines[iline] += '\n'
    for iline in xrange(len(out_w_lines)):
        out_w_lines[iline] += '\n'
    with open(table_name_, 'w') as t:
        t.writelines(out_lines)
        t.write('\n')
        t.writelines(out_w_lines)
        t.close()


def get_ordered_list_of_masses(signal_list_):
 
    tree_mass_total_tmp = []
    tree_mass_total = []
    signal_tree_mass_sorted = OrderedDict()

    for signal in signal_list_:
        signal_tree_mass_sorted[signal] = OrderedDict()

        tree_name_mass = [(int(mass.split('_')[1]), int(mass.split('_')[2])) for mass in signal_list_[signal]['trees']]
        tree_mass_total_tmp.append(tree_name_mass)
        tree_name_mass.sort(key=lambda x: int(x[0]))

        for tree in tree_name_mass:
            if str(tree[0]) not in signal_tree_mass_sorted[signal].keys():
                signal_tree_mass_sorted[signal][str(tree[0])] = []

        for tree in tree_name_mass:
            signal_tree_mass_sorted[signal][str(tree[0])].append(tree[1])

        for stop in signal_tree_mass_sorted[signal]:
            signal_tree_mass_sorted[signal][stop].sort(key=lambda x: int(x))        

    for signal in signal_tree_mass_sorted:
        print signal

        for stop in signal_tree_mass_sorted[signal]:
            print '   ', stop

            for lsp in signal_tree_mass_sorted[signal][stop]:
                print '       ', lsp

    for tree_list in tree_mass_total_tmp:
        for tree in tree_list:
            tree_mass_total.append(tree)

    tree_mass_total.sort(key=lambda x: int(x[0]))
    tree_mass_sorted = OrderedDict()

    for tree in tree_mass_total:
        if tree[0] not in tree_mass_sorted.keys():
            tree_mass_sorted[str(tree[0])] = []

    for tree in tree_mass_total:
        tree_mass_sorted[str(tree[0])].append(tree[1])

    for stop in tree_mass_sorted:
        tree_mass_sorted[stop].sort(key=lambda x: int(x))
       
    return tree_mass_sorted

def make_table_of_masses(ordered_list_):
    """
       input
           @ordered_list_ : input table produced by function 'get_ordered_list_of_masses()'
       output
           @nice_table : 2d table of mass points that exist
    """
     
    lsp_list = []
    nice_table = []
    nice_table.append('{:^7}'.format(' '))

    for stop in ordered_list_:
        for lsp in ordered_list_[stop]:
            if lsp not in lsp_list:
                lsp_list.append(lsp)
        nice_table.append('{:^7}'.format(stop))

    for lsp in lsp_list:
        nice_table[0] += '{:^7}'.format(str(lsp))

    existence_table = np.zeros((len(ordered_list_), len(lsp_list)), dtype=int)

    for istop, stop in enumerate(ordered_list_):
        for ilsp, lsp in enumerate(lsp_list):
            if lsp in ordered_list_[stop]:
                existence_table[istop][ilsp] += 1

    for irow, row in enumerate(existence_table):
        for entry in row:
            nice_table[irow + 1] += '{:^7}'.format(entry)

    for irow in xrange(len(nice_table)):
        nice_table[irow] += '\n'

    with open('table_of_masses.txt', 'w') as t:
        t.writelines(nice_table)
        t.close()
        
     


