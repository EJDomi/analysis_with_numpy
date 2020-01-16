#!/usr/bin/env python

"""
New thing to try out, doing analysis with numpy, converting back to ROOT to do histogramming
Creator: Erich Schmitz
Date: Feb 22, 2019
"""

from time import strftime, localtime
import argparse as arg
import ROOT as rt
import numpy as np
import root_numpy as rnp
import numpy.lib.recfunctions as rfc
import os
from file_table_functions import *
from collections import OrderedDict
import time
import pickle

rt.gROOT.SetBatch()
rt.TH1.AddDirectory(rt.kFALSE)

date = strftime('%d%b%y', localtime())
########## get_histograms function template ############

########################################################


def get_histograms(list_of_files_, variable_list_, cuts_to_apply_=None):
    
    hist = OrderedDict()
    counts = OrderedDict()
    for sample in list_of_files_:
        hist[sample] = OrderedDict()
        counts[sample] = OrderedDict()
        for tree_name in list_of_files_[sample]['trees']:
            print '\nReserving Histograms for:', sample, tree_name 
            hist[sample][tree_name] = OrderedDict()
            counts[sample][tree_name] = OrderedDict()
            # Reserve histograms

            selections_str = [
                 'risr_l_0p5',
                 'risr_0p5',
                 'risr_0p6',
                 'risr_0p7',
                 'risr_0p8',
                 'risr_0p9',
                 'risr_0p95',
                 #'ssv_risr_l_0p5',
                 #'ssv_risr_0p5',
                 #'ssv_risr_0p6',
                 #'ssv_risr_0p7',
                 #'ssv_risr_0p8',
                 #'ssv_risr_0p9',
                 'ssv_risr_0p95',
                 'sbjet_risr_0p95',
                 'stot_risr_0p95',
                ]

            for sel in selections_str:
                hist[sample][tree_name]['PS_'+sel] = rt.TH1D('PS_'+sel+'_'+sample+'_'+tree_name, 'P S system', 512, 0, 1000)
                hist[sample][tree_name]['PT_lep_'+sel] = rt.TH1D('PT_lep_'+sel+'_'+sample+'_'+tree_name, 'PT leps', 512, 0, 1000)
                hist[sample][tree_name]['dphi_lep_S_'+sel] = rt.TH1D('dphi_lep_S_'+sel+'_'+sample+'_'+tree_name, 'PT leps', 512, 0, np.pi)
                hist[sample][tree_name]['dphi_lep_MET_'+sel] = rt.TH1D('dphi_lep_MET_'+sel+'_'+sample+'_'+tree_name, 'PT leps', 512, 0, np.pi)
                hist[sample][tree_name]['PT_lep_dphi_'+sel] = rt.TH2D('PT_lep_dphi_'+sel+'_'+sample+'_'+tree_name, 'PT leps', 512, 0, 1000, 512, 0, np.pi)
                hist[sample][tree_name]['PT_lep_dphi_MET_'+sel] = rt.TH2D('PT_lep_dphi_MET_'+sel+'_'+sample+'_'+tree_name, 'PT leps', 512, 0, 1000, 512, 0, np.pi)
                hist[sample][tree_name]['PS_dphi_'+sel] = rt.TH2D('PS_dphi_'+sel+'_'+sample+'_'+tree_name, 'PT leps', 512, 0, 1000, 512, 0, np.pi)
                hist[sample][tree_name]['PS_dphi_MET_'+sel] = rt.TH2D('PS_dphi_MET_'+sel+'_'+sample+'_'+tree_name, 'PT leps', 512, 0, 1000, 512, 0, np.pi)
                hist[sample][tree_name]['PT_lep_PS_'+sel] = rt.TH2D('PT_lep_PS_'+sel+'_'+sample+'_'+tree_name, 'PT leps', 512, 0, 1000, 512, 0, 1000)
        
    
        for ifile, in_file in enumerate(list_of_files_[sample]['files']):
            for tree_name in list_of_files_[sample]['trees']:
                sample_array = get_tree_info_singular(sample, in_file, tree_name, variable_list_, cuts_to_apply_)
                if sample_array is None: continue

                print '\nGetting Histograms for:', sample, tree_name, in_file
                print 'file: ', ifile+1, ' / ', len(list_of_files_[sample]['files'])

                met = np.array(sample_array['MET'])
                met_phi = np.array(sample_array['MET_phi'])
                ps = np.array(sample_array['PS'])
                pv = np.array(sample_array['PV'])
                risr = np.array(sample_array['RISR'])
                if len(risr) == 0: continue 
                ptisr = np.array(sample_array['PTISR'])

                dphi_lep_s = np.array(sample_array['dphi_lep_S'])

                nlep = np.array(sample_array['Nlep'])
                phi_lep = np.array(sample_array['Phi_lep'])
                nsv = np.array(sample_array['NSV'])
                nsv_s = np.array(sample_array['NSV_S'])
                nsv_isr = np.array(sample_array['NSV_ISR'])
                nbjet_s = np.array(sample_array['Nbjet_S'])
                nbjet = np.array(sample_array['Nbjet'])
                nbjet_isr = np.array(sample_array['Nbjet_ISR'])
                eta_sv = np.array(sample_array['Eta_SV'])
                isr_index_sv = np.array(sample_array['index_SV_ISR'])
                s_index_sv = np.array(sample_array['index_SV_S'])

                weight = np.array(sample_array['weight'])

                id_lep = np.array(sample_array['ID_lep'])
                pt_lep = np.array(sample_array['PT_lep'])
                mini_lep = np.array(sample_array['MiniIso_lep'])
                sip3d_lep = np.array(sample_array['SIP3D_lep'])
                
                risr = np.array([entry[:2] for entry in risr])
                ps = np.array([entry[:2] for entry in ps])
                pv = np.array([entry[:2] for entry in pv])
                ptisr = np.array([entry[:2] for entry in ptisr])
                dphi_lep_s = np.array([entry[:2] for entry in dphi_lep_s])

                nsv_s = np.array([entry[:2] for entry in nsv_s])
                nsv_isr = np.array([entry[:2] for entry in nsv_isr])
                nbjet_s = np.array([entry[:2] for entry in nbjet_s])
                nbjet_isr = np.array([entry[:2] for entry in nbjet_isr])
                isr_index_sv = np.array([entry[:2] if np.shape(entry[:2]) != (0,) else [np.array([],dtype=np.int32),np.array([],dtype=np.int32)] for entry in isr_index_sv])
                s_index_sv = np.array([entry[:2] if np.shape(entry[:2]) != (0,) else [np.array([],dtype=np.int32),np.array([],dtype=np.int32)] for entry in s_index_sv])

                risr = risr[:, 1]
                ps = ps[:, 1]
                pv = pv[:, 1]
                ptisr = ptisr[:, 1]
                dphi_lep_s = dphi_lep_s[:, 1]
                isr_index_sv = isr_index_sv[:, 1]
                s_index_sv = s_index_sv[:, 1]


                nsv_s = nsv_s[:, 1]
                nsv_isr = nsv_isr[:, 1]
                nbjet_s = nbjet_s[:, 1]
                nbjet_isr = nbjet_isr[:, 1]

                weight = 137. * weight
               
                eta_sv_s = np.array([eta[ind] for eta, ind in zip(eta_sv, s_index_sv)])
                eta_sv_isr = np.array([eta[ind] for eta, ind in zip(eta_sv, isr_index_sv)])
               
                nsv = np.array([len(eta[np.abs(eta)<1.5]) for eta in eta_sv])
                nsv_s = np.array([len(eta[np.abs(eta)<1.5]) for eta in eta_sv_s])
                nsv_isr = np.array([len(eta[np.abs(eta)<1.5]) for eta in eta_sv_isr])

 
                one_lep = nlep == 1 
                two_lep = nlep == 2 
                ################        Medium       #####################
                mini_lep = np.array([mini[lid>=3]*pt[lid>=3] for mini, lid, pt in zip(mini_lep, id_lep, pt_lep)])
                sip3d_lep = np.array([ips[lid>=3] for ips, lid in zip(sip3d_lep, id_lep)])
                phi_lep = np.array([ips[lid>=3] for ips, lid in zip(phi_lep, id_lep)])

                #pt_med_lep = np.array([lep[lid>=2] for lep, lid in zip(pt_lep, id_lep)]) 
                #pt_mini_lep = np.array([lep[mini<6] for lep, mini in zip(pt_med_lep, mini_lep)])
                #eta_med_lep = np.array([lep[lid>=2] for lep, lid in zip(eta_lep, id_lep)]) 
                #eta_mini_lep = np.array([lep[mini<6] for lep, mini in zip(eta_med_lep, mini_lep)])
                #phi_med_lep = np.array([lep[lid>=2] for lep, lid in zip(phi_lep, id_lep)]) 
                #phi_mini_lep = np.array([lep[mini<6] for lep, mini in zip(phi_med_lep, mini_lep)])
                #m_med_lep = np.array([lep[lid>=2] for lep, lid in zip(m_lep, id_lep)]) 
                #m_mini_lep = np.array([lep[mini<0.1] for lep, mini in zip(m_med_lep, mini_lep)])

                phi_mini_lep = np.array([ phi[np.logical_and(sip<4., mini<6.0)] for phi, sip, mini in zip(phi_lep, sip3d_lep, mini_lep)])
                dphi_lep_met = np.array([np.array([phi - m_phi for phi in leps])
                                            for leps, m_phi in zip(phi_mini_lep, met_phi)])

                dphi_lep_met = np.array([ np.array([phi + 2*np.pi if phi < -np.pi else phi for phi in leps]) for leps in dphi_lep_met])
                dphi_lep_met = np.array([ np.array([phi - 2*np.pi if phi >= np.pi else phi for phi in leps]) for leps in dphi_lep_met])
                dphi_lep_met = np.abs(dphi_lep_met) 

                mini_med_mask = np.array([ True if len(mini[np.logical_and(mini<6, sip<4.)]) == 1 and one else False for mini, sip, one in zip(mini_lep, sip3d_lep, one_lep)])
                #mini_med_mask = np.array([ True if len(mini[np.logical_and(mini<6, sip3d<4.)]) == 1 and one else False for mini, sip3d, one in zip(mini_lep,sip3d_lep,one_lep)])
                #mini_med_mask = np.array([ True if len(mini[np.logical_and(mini<6, sip3d<4.)]) == 2 and two else False for mini, sip3d, two in zip(mini_lep,sip3d_lep,two_lep)])
               
 
                risr_0p95 = risr > 0.95
                risr_0p9 = risr > 0.9
                risr_0p8 = risr > 0.8
                risr_0p7 = risr > 0.7
                risr_0p6 = risr > 0.6
                risr_0p5 = risr > 0.5

                risr_l_0p5 = risr <= 0.5
                risr_l_0p6 = risr <= 0.6
                risr_l_0p7 = risr <= 0.7
                risr_l_0p8 = risr <= 0.8
                risr_l_0p9 = risr <= 0.9
                risr_l_0p95 = risr <= 0.95

                met_200 = met > 200
                ptisr_200 = ptisr > 200
                no_isrb = nbjet_isr < 1
                ssv_ge1 = nsv_s > 0
                sbjet_ge1 = nbjet_s > 0
                stot_ge1 = nsv_s + nbjet_s > 0

                evt_ssv_risr_l_0p5_selection_mask = np.all([met_200, risr_l_0p5, ptisr_200, mini_med_mask, no_isrb, ssv_ge1], axis=0)
                evt_ssv_risr_0p5_selection_mask = np.all([met_200, risr_0p5, risr_l_0p6, ptisr_200, mini_med_mask, no_isrb, ssv_ge1], axis=0)
                evt_ssv_risr_0p6_selection_mask = np.all([met_200, risr_0p6, risr_l_0p7, ptisr_200, mini_med_mask, no_isrb, ssv_ge1], axis=0)
                evt_ssv_risr_0p7_selection_mask = np.all([met_200, risr_0p7, risr_l_0p8, ptisr_200, mini_med_mask, no_isrb, ssv_ge1], axis=0)
                evt_ssv_risr_0p8_selection_mask = np.all([met_200, risr_0p8, risr_l_0p9, ptisr_200, mini_med_mask, no_isrb, ssv_ge1], axis=0)
                evt_ssv_risr_0p9_selection_mask = np.all([met_200, risr_0p9, risr_l_0p95, ptisr_200, mini_med_mask, no_isrb, ssv_ge1], axis=0)
                evt_ssv_risr_0p95_selection_mask = np.all([met_200, risr_0p95, ptisr_200, mini_med_mask, no_isrb, ssv_ge1], axis=0)

                evt_risr_l_0p5_selection_mask = np.all([met_200, risr_l_0p5, ptisr_200, mini_med_mask, no_isrb], axis=0)
                evt_risr_0p5_selection_mask = np.all([met_200, risr_0p5, risr_l_0p6, ptisr_200, mini_med_mask, no_isrb], axis=0)
                evt_risr_0p6_selection_mask = np.all([met_200, risr_0p6, risr_l_0p7, ptisr_200, mini_med_mask, no_isrb], axis=0)
                evt_risr_0p7_selection_mask = np.all([met_200, risr_0p7, risr_l_0p8, ptisr_200, mini_med_mask, no_isrb], axis=0)
                evt_risr_0p8_selection_mask = np.all([met_200, risr_0p8, risr_l_0p9, ptisr_200, mini_med_mask, no_isrb], axis=0)
                evt_risr_0p9_selection_mask = np.all([met_200, risr_0p9, risr_l_0p95, ptisr_200, mini_med_mask, no_isrb], axis=0)
                evt_risr_0p95_selection_mask = np.all([met_200, risr_0p95, ptisr_200, mini_med_mask, no_isrb], axis=0)

                evt_sbjet_risr_0p95_selection_mask = np.all([met_200, risr_0p95, ptisr_200, mini_med_mask, no_isrb, sbjet_ge1], axis=0)
                evt_stot_risr_0p95_selection_mask = np.all([met_200, risr_0p95, ptisr_200, mini_med_mask, no_isrb, stot_ge1], axis=0)


                selections_str = [
                                  'risr_l_0p5',
                                  'risr_0p5',
                                  'risr_0p6',
                                  'risr_0p7',
                                  'risr_0p8',
                                  'risr_0p9',
                                  'risr_0p95',
                                  #'ssv_risr_l_0p5',
                                  #'ssv_risr_0p5',
                                  #'ssv_risr_0p6',
                                  #'ssv_risr_0p7',
                                  #'ssv_risr_0p8',
                                  #'ssv_risr_0p9',
                                  'ssv_risr_0p95',
                                  'sbjet_risr_0p95',
                                  'stot_risr_0p95',
                                 ]

                selections = [
                              evt_risr_l_0p5_selection_mask,
                              evt_risr_0p5_selection_mask,
                              evt_risr_0p6_selection_mask,
                              evt_risr_0p7_selection_mask,
                              evt_risr_0p8_selection_mask,
                              evt_risr_0p9_selection_mask,
                              evt_risr_0p95_selection_mask,
                              #evt_ssv_risr_l_0p5_selection_mask,
                              #evt_ssv_risr_0p5_selection_mask,
                              #evt_ssv_risr_0p6_selection_mask,
                              #evt_ssv_risr_0p7_selection_mask,
                              #evt_ssv_risr_0p8_selection_mask,
                              #evt_ssv_risr_0p9_selection_mask,
                              evt_ssv_risr_0p95_selection_mask,
                              evt_sbjet_risr_0p95_selection_mask,
                              evt_stot_risr_0p95_selection_mask,
                             ]
          
                for sel, sel_str in zip(selections, selections_str): 
                    if np.any(sel):

                        tmp_weight = weight[sel]

                        pt_lep_sel = np.array([pt for pt, select in zip(pt_lep, sel) if select])
                        dphi_lep_met_sel = np.array([pt for pt, select in zip(dphi_lep_met, sel) if select])
                        #tmp_weight_lep = np.array([ [w]*len(lep) for w, lep in zip(tmp_weight, pt_lep_sel)])

                        ps_sel = ps[sel]

                        dphi_lep_s_sel = dphi_lep_s[sel]

                        pt_lep_sel = np.concatenate(pt_lep_sel)
                        dphi_lep_met_sel = np.concatenate(dphi_lep_met_sel)
                        #tmp_weight_lepl = np.concatenate(tmp_weight_lep)

                        rnp.fill_hist(hist[sample][tree_name]['PS_'+sel_str], ps_sel, tmp_weight) 
                        rnp.fill_hist(hist[sample][tree_name]['PT_lep_'+sel_str], pt_lep_sel, tmp_weight) 
                        rnp.fill_hist(hist[sample][tree_name]['dphi_lep_S_'+sel_str], dphi_lep_s_sel, tmp_weight) 
                        rnp.fill_hist(hist[sample][tree_name]['dphi_lep_MET_'+sel_str], dphi_lep_met_sel, tmp_weight) 
                        rnp.fill_hist(hist[sample][tree_name]['PT_lep_dphi_'+sel_str], np.swapaxes([pt_lep_sel, dphi_lep_s_sel],0,1), tmp_weight) 
                        rnp.fill_hist(hist[sample][tree_name]['PT_lep_dphi_MET_'+sel_str], np.swapaxes([pt_lep_sel, dphi_lep_met_sel],0,1), tmp_weight) 
                        rnp.fill_hist(hist[sample][tree_name]['PS_dphi_'+sel_str], np.swapaxes([ps_sel, dphi_lep_s_sel],0,1), tmp_weight) 
                        rnp.fill_hist(hist[sample][tree_name]['PS_dphi_MET_'+sel_str], np.swapaxes([ps_sel, dphi_lep_met_sel],0,1), tmp_weight) 
                        rnp.fill_hist(hist[sample][tree_name]['PT_lep_PS_'+sel_str], np.swapaxes([pt_lep_sel, ps_sel],0,1), tmp_weight) 
                    print 'finished filling: ' + sel_str
                print 'finished filling'
    return hist


             

if __name__ == "__main__":

    parser = arg.ArgumentParser(description='receiving file lists for batch jobs')
    parser.add_argument('-f', '--file_list', type=str)
    parser.add_argument('-o', '--out_file', type=str) 

    args = parser.parse_args()

    sample_file = args.file_list
    file_name = args.out_file

    sample_list = pickle.load( open(sample_file, "rb"))
    variables = ['MET', 'MET_phi', 'RISR', 'PTISR', 'PT_lep', 'Phi_lep', 'MiniIso_lep', 'ID_lep', 'Nbjet', 'Nbjet_S', 'Nbjet_ISR', 'NSV_S', 'NSV_ISR', 'weight', 'NSV', 'Nlep', 'Eta_SV', 'index_SV_S','index_SV_ISR', 'SIP3D_lep', 'PS', 'PV', 'PLa', 'PLb', 'dphi_lep_S']

    start_b = time.time()    

    hist_sample = get_histograms(sample_list, variables, None)

    write_hists_to_file(hist_sample, file_name) 
    stop_b = time.time()
    print "total:      ", stop_b - start_b


