import ROOT as rt
import root_numpy as rnp
from time import strftime, localtime
from collections import OrderedDict
from plotting_susy_cff import plot_configurables as pc
from plotting_susy_cff import sample_configurables as sc
import imp, os
from file_table_functions import *

date = strftime('%d%b%y', localtime())


def read_in_hists(in_file_):
    print 'look at all that posterity' 
    in_file = rt.TFile(in_file_, 'r')
    hists = OrderedDict()
    for key in in_file.GetListOfKeys():
        print key.GetName()
        key_name = key.GetName()
        if key_name not in sc: continue
        sample_dir = in_file.Get(key.GetName())
        hists[key_name] = OrderedDict()
        for tree_key in sample_dir.GetListOfKeys():
            print tree_key.GetName()
            tree_name = tree_key.GetName()
            tree_dir = sample_dir.Get(tree_key.GetName())
            hists[key_name][tree_name] = OrderedDict()
            for hist_key in tree_dir.GetListOfKeys():
                print hist_key.GetName()
                hist_name = hist_key.GetName()
                hist = tree_dir.Get(hist_key.GetName())
                if 'DY' in key_name:
                    hists[key_name][tree_name]['_'.join(hist_name.split('_')[:-4])] = hist
                elif 'SMS_T2bW' in key_name:
                    hists[key_name][tree_name]['_'.join(hist_name.split('_')[:-4])] = hist
                elif 'SMS_T2_' in key_name:
                    hists[key_name][tree_name]['_'.join(hist_name.split('_')[:-5])] = hist
                elif 'SMS' in key_name:
                    hists[key_name][tree_name]['_'.join(hist_name.split('_')[:-5])] = hist
		elif '2017' in key_name:
                    hists[key_name][tree_name]['_'.join(hist_name.split('_')[:-3])] = hist
                else:
                    hists[key_name][tree_name]['_'.join(hist_name.split('_')[:-2])] = hist
    return hists

def make_nsv_table(sig_hists_, b_hists_):

    pre_pend = 'N,ttjets,wjets,T2-4bd_500_420,T2-4bd_500_490,T2bW_500_460'

    out_lines = OrderedDict()
    out_lines['SVs'] = []
    out_lines['bjets'] = []
    out_lines['Bs'] = []
    
    hist_objects = ['SVs', 'bjets', 'Bs']
    name = 'N_{}'
    name_S = 'NS_{}'
    name_ISR = 'NISR_{}'

    for obj in hist_objects:
        ttjets_inc = rnp.hist2array(b_hists_['TTJets_2017']['KUAnalysis'][name.format(obj)])
        ttjets_s = rnp.hist2array(b_hists_['TTJets_2017']['KUAnalysis'][name_S.format(obj)])
        ttjets_isr = rnp.hist2array(b_hists_['TTJets_2017']['KUAnalysis'][name_ISR.format(obj)])

        wjets_inc = rnp.hist2array(b_hists_['WJets_2017']['KUAnalysis'][name.format(obj)])
        wjets_s = rnp.hist2array(b_hists_['WJets_2017']['KUAnalysis'][name_S.format(obj)])
        wjets_isr = rnp.hist2array(b_hists_['WJets_2017']['KUAnalysis'][name_ISR.format(obj)])

        t420_inc = rnp.hist2array(sig_hists_['SMS-T2-4bd_420']['SMS_500_420'][name.format(obj)])
        t420_s = rnp.hist2array(sig_hists_['SMS-T2-4bd_420']['SMS_500_420'][name_S.format(obj)])
        t420_isr = rnp.hist2array(sig_hists_['SMS-T2-4bd_420']['SMS_500_420'][name_ISR.format(obj)])

        t490_inc = rnp.hist2array(sig_hists_['SMS-T2-4bd_490']['SMS_500_490'][name.format(obj)])
        t490_s = rnp.hist2array(sig_hists_['SMS-T2-4bd_490']['SMS_500_490'][name_S.format(obj)])
        t490_isr = rnp.hist2array(sig_hists_['SMS-T2-4bd_490']['SMS_500_490'][name_ISR.format(obj)])

        t460_inc = rnp.hist2array(sig_hists_['SMS-T2bW_dM']['SMS_500_460'][name.format(obj)])
        t460_s = rnp.hist2array(sig_hists_['SMS-T2bW_dM']['SMS_500_460'][name_S.format(obj)])
        t460_isr = rnp.hist2array(sig_hists_['SMS-T2bW_dM']['SMS_500_460'][name_ISR.format(obj)])

        out_inc = ['{},{},{},{},{},{}'.format(ibin, tt, w, s42, s49, s46) 
                   for ibin, tt, w, s42, s49, s46 in 
                   zip(xrange(len(ttjets_inc)), ttjets_inc, wjets_inc, t420_inc, t490_inc, t460_inc)]
        
        out_s = ['{},{},{},{},{},{}'.format(ibin, tt, w, s42, s49, s46) 
                   for ibin, tt, w, s42, s49, s46 in 
                   zip(xrange(len(ttjets_s)), ttjets_s, wjets_s, t420_s, t490_s, t460_s)]
        
        out_isr = ['{},{},{},{},{},{}'.format(ibin, tt, w, s42, s49, s46) 
                   for ibin, tt, w, s42, s49, s46 in 
                   zip(xrange(len(ttjets_isr)), ttjets_isr, wjets_isr, t420_isr, t490_isr, t460_isr)]

        out_lines[obj].append(pre_pend)
        out_lines[obj].append(name.format(obj))
        for out in out_inc:
            out_lines[obj].append(out)

        out_lines[obj].append(name_S.format(obj))
        for out in out_s:
            out_lines[obj].append(out)

        out_lines[obj].append(name_ISR.format(obj))
        for out in out_isr:
            out_lines[obj].append(out)
        
        for ilines, out in enumerate(out_lines[obj]):
           out_lines[obj][ilines] += '\n'

        with open('table_'+obj+'.csv', 'w') as t:
            t.writelines(out_lines[obj])

if __name__ == "__main__":

    signal_file = './output_signal_risr_0p8_mixed_nsv_19Nov19.root'
    sig_hists = read_in_hists(signal_file)
    background_file = './output_background_risr_0p8_mixed_nsv_18Nov19.root'
    b_hists = read_in_hists(background_file)

    make_nsv_table(sig_hists, b_hists)
