#!/usr/bin/env python

"""
New thing to try out, doing analysis with numpy, converting back to ROOT to do histogramming
Creator: Erich Schmitz
Date: Feb 22, 2019
"""

from time import strftime, localtime
import ROOT as rt
import numpy as np
import root_numpy as rnp
import numpy.lib.recfunctions as rfc
import os
from file_table_functions import *
from collections import OrderedDict
import time
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
            hist[sample][tree_name]['MET'] = rt.TH1D('MET_'+sample+'_'+tree_name, 'E_{T}^{miss} [GeV]', 500, 0, 1000)
            hist[sample][tree_name]['MET_postrisr'] = rt.TH1D('MET_postrisr_'+sample+'_'+tree_name, 'E_{T}^{miss} [GeV]', 1024, 0, 1000)
            hist[sample][tree_name]['NISR_bjets'] = rt.TH1D('NISR_bjets_'+sample+'_'+tree_name, 'N ISR b jets', 500, 0, 20)
            hist[sample][tree_name]['NS_jets'] = rt.TH1D('NS_jets_'+sample+'_'+tree_name, 'N S jets', 500, 0, 20)
            hist[sample][tree_name]['NS_SVs'] = rt.TH1D('NS_SVs_'+sample+'_'+tree_name, 'N S SVs', 500, 0, 20)
            hist[sample][tree_name]['NS_STs'] = rt.TH1D('NS_STs_'+sample+'_'+tree_name, 'N S STs', 500, 0, 20)
            hist[sample][tree_name]['NS_SBs'] = rt.TH1D('NS_SBs_'+sample+'_'+tree_name, 'N S SBs', 500, 0, 20)
            hist[sample][tree_name]['dphi_lep_MET'] = rt.TH1D('dphi_lep_MET_'+sample+'_'+tree_name, 'dphi l met', 1024, -np.pi, np.pi)
            hist[sample][tree_name]['MS'] = rt.TH1D('MS_'+sample+'_'+tree_name, 'MS', 1024, 0, 1000)
            hist[sample][tree_name]['MS_preisrb'] = rt.TH1D('MS_preisrb_'+sample+'_'+tree_name, 'MS pre isrb cut', 1024, 0, 1000)
            hist[sample][tree_name]['MS_prest'] = rt.TH1D('MS_prest_'+sample+'_'+tree_name, 'MS pre st cut', 1024, 0, 1000)
            hist[sample][tree_name]['PzS_div_PTISR'] = rt.TH1D('PzS_div_PTISR_'+sample+'_'+tree_name, 'PzS / PzS + MS', 1024, 0, 1)
            hist[sample][tree_name]['PTS_div_PTISR'] = rt.TH1D('PTS_div_PTISR_'+sample+'_'+tree_name, 'PTS / PTS + MS', 1024, 0, 1)
            hist[sample][tree_name]['RISR_PTISR'] = rt.TH2D('RISR_PTISR_'+sample+'_'+tree_name, 'MS vs RISR', 1024, 0, 2, 1024, 0, 1000)
            hist[sample][tree_name]['RISR_PTISR_postscut'] = rt.TH2D('RISR_PTISR_postscut_'+sample+'_'+tree_name, 'MS vs RISR', 1024, 0, 2, 1024, 0, 1000)
            hist[sample][tree_name]['RISR_PTISR_prest'] = rt.TH2D('RISR_PTISR_prest_'+sample+'_'+tree_name, 'MS vs RISR', 1024, 0, 2, 1024, 0, 1000)
            hist[sample][tree_name]['RISR_PTISR_preisrb'] = rt.TH2D('RISR_PTISR_preisrb_'+sample+'_'+tree_name, 'MS vs RISR', 1024, 0, 2, 1024, 0, 1000)
           
            hist[sample][tree_name]['cos_lep_S'] = rt.TH1D('cos_lep_S_'+sample+'_'+tree_name, 'cos_lep_S', 1024, -1, 1)
            hist[sample][tree_name]['dphi_lep_S'] = rt.TH1D('dphi_lep_S_'+sample+'_'+tree_name, 'dphi_lep_S', 1024, -np.pi, np.pi)
            
        for ifile, in_file in enumerate(list_of_files_[sample]['files']):
            for tree_name in list_of_files_[sample]['trees']:
                sample_array = get_tree_info_singular(sample, in_file, tree_name, variable_list_, cuts_to_apply_)
                if sample_array is None: continue

                print '\nGetting Histograms for:', sample, tree_name, in_file
                print 'file: ', ifile+1, ' / ', len(list_of_files_[sample]['files'])

                met = np.array(sample_array['MET'])
                phi_met = np.array(sample_array['MET_phi'])
                risr = np.array(sample_array['RISR'])
                ptisr = np.array(sample_array['PTISR'])

                ms = np.array(sample_array['MS'])
                nsv_s = np.array(sample_array['NSV_S'])
                nbjet_s = np.array(sample_array['Nbjet_S'])
                njet_s = np.array(sample_array['Njet_S'])
                nbjet_isr = np.array(sample_array['Nbjet_ISR'])
                ps = np.array(sample_array['PS'])
                pts = np.array(sample_array['PTS'])
                pzs = np.array(sample_array['PzS'])
                cos_lep_s = np.array(sample_array['cos_lep_S'])
                dphi_lep_s = np.array(sample_array['dphi_lep_S'])

                weight = np.array(sample_array['weight'])

                id_lep = np.array(sample_array['ID_lep'])
                pt_lep = np.array(sample_array['PT_lep'])
                phi_lep = np.array(sample_array['Phi_lep'])
                eta_lep = np.array(sample_array['Eta_lep'])
                m_lep = np.array(sample_array['M_lep'])
                mini_lep = np.array(sample_array['MiniIso_lep'])
                
                risr = np.array([entry[:2] for entry in risr])
                ptisr = np.array([entry[:2] for entry in ptisr])

                ms = np.array([entry[:2] for entry in ms])
                nsv_s = np.array([entry[:2] for entry in nsv_s])
                nbjet_s = np.array([entry[:2] for entry in nbjet_s])
                njet_s = np.array([entry[:2] for entry in njet_s])
                nbjet_isr = np.array([entry[:2] for entry in nbjet_isr])
                ps = np.array([entry[:2] for entry in ps])
                pts = np.array([entry[:2] for entry in pts])
                pzs = np.array([entry[:2] for entry in pzs])
                cos_lep_s = np.array([entry[:2] for entry in cos_lep_s])
                dphi_lep_s = np.array([entry[:2] for entry in dphi_lep_s])

                risr = risr[:, 1]
                ptisr = ptisr[:, 1]

                ms = ms[:, 1]
                nsv_s = nsv_s[:, 1]
                nbjet_s = nbjet_s[:, 1]
                njet_s = njet_s[:, 1]
                nbjet_isr = nbjet_isr[:, 1]
                ps = ps[:, 1]
                pts = pts[:, 1]
                pzs = pzs[:, 1]
                cos_lep_s = cos_lep_s[:, 1]
                dphi_lep_s = dphi_lep_s[:, 1]

                weight = 137. * weight
                

                ################        Medium       #####################
                mini_lep = np.array([mini[lid>=2] for mini, lid in zip(mini_lep, id_lep)])

                pt_med_lep = np.array([lep[lid>=2] for lep, lid in zip(pt_lep, id_lep)]) 
                pt_mini_lep = np.array([lep[mini<0.1] for lep, mini in zip(pt_med_lep, mini_lep)])
                eta_med_lep = np.array([lep[lid>=2] for lep, lid in zip(eta_lep, id_lep)]) 
                eta_mini_lep = np.array([lep[mini<0.1] for lep, mini in zip(eta_med_lep, mini_lep)])
                phi_med_lep = np.array([lep[lid>=2] for lep, lid in zip(phi_lep, id_lep)]) 
                phi_mini_lep = np.array([lep[mini<0.1] for lep, mini in zip(phi_med_lep, mini_lep)])
                m_med_lep = np.array([lep[lid>=2] for lep, lid in zip(m_lep, id_lep)]) 
                m_mini_lep = np.array([lep[mini<0.1] for lep, mini in zip(m_med_lep, mini_lep)])

                
                risr_0p95 = risr > 0.95
                risr_0p8 = risr > 0.8
                met_200 = met > 200
                ptisr_200 = ptisr > 200
                ms_80 = ms < 80
                ms_120 = ms > 120

                sj_ge_1 = njet_s >= 1
                st_ge_1 = (nsv_s + njet_s) >= 1
                st_0 = (nsv_s + njet_s) < 1
                sv_ge_1 = nsv_s >= 1
                sj_0 = njet_s < 1
                b_isr_ge_1 = nbjet_isr >= 1
                b_isr_ge_2 = nbjet_isr >= 2
                b_isr_none = nbjet_isr < 1
################## ttbar 1l #########################                
#                mini_med_mask = np.array([ True if len(mini[mini<0.1]) == 1 and len(lid) == 1 else False for mini, lid in zip(mini_lep,id_lep)])
#                evt_risr_selection_mask = np.all([met_200], axis=0)
#                evt_isrb_selection_mask = np.all([met_200, risr_0p8, ptisr_200, mini_med_mask], axis=0)
#                evt_sj_selection_mask = np.all([met_200, risr_0p8, ptisr_200, mini_med_mask, b_isr_ge_2], axis=0)
#                evt_selection_mask = np.all([met_200, risr_0p8, ptisr_200, mini_med_mask, b_isr_ge_2, sv_ge_1], axis=0)
################## ttbar 2l ##########################
#                mini_med_mask = np.array([ True if len(mini[mini<0.1]) == 2 and len(lid) == 2 else False for mini, lid in zip(mini_lep,id_lep)])
#                evt_risr_selection_mask = np.all([met_200], axis=0)
#                evt_isrb_selection_mask = np.all([met_200, risr_0p8, ptisr_200, mini_med_mask], axis=0)
#                evt_sj_selection_mask = np.all([met_200, risr_0p8, ptisr_200, mini_med_mask, b_isr_ge_1], axis=0)
#                evt_selection_mask = np.all([met_200, risr_0p8, ptisr_200, mini_med_mask, b_isr_ge_1, sv_ge_1], axis=0)
################## WJet 1l ###########################
                mini_med_mask = np.array([ True if len(mini[mini<0.1]) == 1 and len(lid) == 1 else False for mini, lid in zip(mini_lep,id_lep)])
                evt_risr_selection_mask = np.all([met_200], axis=0)
                evt_isrb_selection_mask = np.all([met_200, risr_0p8, ptisr_200, mini_med_mask], axis=0)
                evt_sj_selection_mask = np.all([met_200, risr_0p8, ptisr_200, mini_med_mask, b_isr_none], axis=0)
                evt_selection_mask = np.all([met_200, risr_0p8, ptisr_200, mini_med_mask, b_isr_none, st_0], axis=0)

                phi_evt_lep = np.array([lep[0] for lep, mask in zip(phi_mini_lep, evt_selection_mask) if mask])

                met_weight = weight
                risr_weight = weight[evt_risr_selection_mask]
                isrb_weight = weight[evt_isrb_selection_mask]
                sj_weight = weight[evt_sj_selection_mask]
                evt_weight = weight[evt_selection_mask]

                risr_isrb = risr[evt_isrb_selection_mask]
                met_isrb = met[evt_isrb_selection_mask]
                ptisr_risr_isrb = ptisr[evt_isrb_selection_mask]

                risr_sj = risr[evt_sj_selection_mask]
                ptisr_risr_sj = ptisr[evt_sj_selection_mask]

                risr_evt = risr[evt_selection_mask]
                ptisr_risr_evt = ptisr[evt_selection_mask]

                risr = risr[evt_risr_selection_mask]
                ptisr_risr = ptisr[evt_risr_selection_mask]

                nbjet_isr = nbjet_isr[evt_isrb_selection_mask]

                njet_s = njet_s[evt_sj_selection_mask]
                nsv_s = nsv_s[evt_sj_selection_mask]
                nst_s = njet_s + nsv_s
                nsb_s = nbjet_s[evt_sj_selection_mask] + nsv_s

                phi_met = phi_met[evt_selection_mask]

                pts = pts[evt_selection_mask]
                pzs = pzs[evt_selection_mask]
                ptisr = ptisr[evt_selection_mask]
                cos_lep_s = cos_lep_s[evt_selection_mask]
                dphi_lep_s = dphi_lep_s[evt_selection_mask]

                
                ms_isrb = ms[evt_isrb_selection_mask]
                ms_st = ms[evt_sj_selection_mask]
                ms = ms[evt_selection_mask]

##################  delta phi ###########################
                dphi_l_met = phi_evt_lep - phi_met

                dphi_l_met = np.array([ phi + 2*np.pi if phi < -np.pi else phi for phi in dphi_l_met])
                dphi_l_met = np.array([ phi - 2*np.pi if phi >= np.pi else phi for phi in dphi_l_met])
#########################################################

                div_pz_ptisr = np.abs(pzs) / (np.abs(pzs) + ptisr)
                div_pts_ptisr = pts / (pts + ptisr)
                pzs = np.abs(pzs) 
                ps = np.abs(ps) 

                rnp.fill_hist(hist[sample][tree_name]['MET'], met, weight) 
                if np.any(evt_risr_selection_mask):
                    rnp.fill_hist(hist[sample][tree_name]['RISR_PTISR'], np.swapaxes([risr, ptisr_risr],0,1), risr_weight) 

                if np.any(evt_isrb_selection_mask):
                    rnp.fill_hist(hist[sample][tree_name]['RISR_PTISR_preisrb'], np.swapaxes([risr_isrb, ptisr_risr_isrb],0,1), isrb_weight) 
                    rnp.fill_hist(hist[sample][tree_name]['NISR_bjets'], nbjet_isr, isrb_weight) 
                    rnp.fill_hist(hist[sample][tree_name]['MS_preisrb'], ms_isrb, isrb_weight) 
                    rnp.fill_hist(hist[sample][tree_name]['MET_postrisr'], met_isrb, isrb_weight) 

                if np.any(evt_sj_selection_mask):
                    rnp.fill_hist(hist[sample][tree_name]['RISR_PTISR_prest'], np.swapaxes([risr_sj, ptisr_risr_sj],0,1), sj_weight) 
                    rnp.fill_hist(hist[sample][tree_name]['NS_jets'], njet_s, sj_weight) 
                    rnp.fill_hist(hist[sample][tree_name]['NS_SVs'], nsv_s, sj_weight) 
                    rnp.fill_hist(hist[sample][tree_name]['NS_STs'], nst_s, sj_weight) 
                    rnp.fill_hist(hist[sample][tree_name]['NS_SBs'], nsb_s, sj_weight) 
                    rnp.fill_hist(hist[sample][tree_name]['MS_prest'], ms_st, sj_weight) 
                                    
                if not np.any(evt_selection_mask):
                    print 'finished filling'
                    continue
                rnp.fill_hist(hist[sample][tree_name]['RISR_PTISR_postscut'], np.swapaxes([risr_evt, ptisr_risr_evt],0,1), evt_weight) 
                rnp.fill_hist(hist[sample][tree_name]['dphi_lep_MET'], dphi_l_met, evt_weight) 
                rnp.fill_hist(hist[sample][tree_name]['cos_lep_S'], cos_lep_s, evt_weight) 
                rnp.fill_hist(hist[sample][tree_name]['dphi_lep_S'], dphi_lep_s, evt_weight) 
                rnp.fill_hist(hist[sample][tree_name]['MS'], ms, evt_weight) 
                rnp.fill_hist(hist[sample][tree_name]['PzS_div_PTISR'], div_pz_ptisr, evt_weight) 
                rnp.fill_hist(hist[sample][tree_name]['PTS_div_PTISR'], div_pts_ptisr, evt_weight) 

                print 'finished filling'
    return hist


             

if __name__ == "__main__":

    signals = { 
#    'SMS-T2bW_dM' : [
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_SMS_Stop/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_94X',
#], 
    'SMS-T2-4bd_420' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_SMS_Stop/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
#                    '/home/t3-ku/erichjs/work/Ewkinos/reducer/CMSSW_10_1_4_patch1/src/KUEWKinoAnalysis_dev_v2/output_10Sep19/Fall17_94X_SMS_I/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
],
    'SMS-T2-4bd_490' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_SMS_Stop/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
#                    '/home/t3-ku/erichjs/work/Ewkinos/reducer/CMSSW_10_1_4_patch1/src/KUEWKinoAnalysis_dev_v2/output_10Sep19/Fall17_94X_SMS_I/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
],
#    'SMS-TChiWH' : [ 
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_SMS/SMS-TChiWH_WToLNu_HToBB_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_94X',
#],
              }
    backgrounds = {
    'TTJets_2017' : [
#                     '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X'
#                     '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X'
                     '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkgextra/TTTT_TuneCP5_PSweights_13TeV-amcatnlo-pythia8_Fall17_94X',
                     '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkgextra/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                     '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkgextra/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                     '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkgextra/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
],
    'WJets_2017' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
],
    'ZJets_2017' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-100To200_13TeV-madgraph_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-200To400_13TeV-madgraph_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-400To600_13TeV-madgraph_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-600To800_13TeV-madgraph_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-800To1200_13TeV-madgraph_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph_Fall17_94X',
],
    'DY_M50_2017' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
],
    'WW_2017' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WWTo4Q_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WWG_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_94X',
],
    'ZZ_2017' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZZTo2L2Nu_13TeV_powheg_pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZZTo2Q2Nu_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZZTo4L_13TeV_powheg_pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_94X',
],
    'WZ_2017' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZG_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZZ_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZ_TuneCP5_13TeV-pythia8_Fall17_94X',
],
                  }
    variables = ['MET', 'MET_phi', 'MS', 'PS', 'PTS', 'PzS', 'dphiSI', 'dphiCMI', 'RISR', 'PTISR', 'PT_lep', 'Eta_lep', 'Phi_lep', 'M_lep', 'MiniIso_lep', 'ID_lep', 'Njet_S', 'Nbjet_S', 'Nbjet_ISR', 'NSV_S', 'cos_lep_S', 'dphi_lep_S', 'weight']

    start_b = time.time()    
    background_list = process_the_samples(backgrounds, None, None)
    hist_background = get_histograms(background_list, variables, None)

    write_hists_to_file(hist_background, './output_background_wjet_1l_control_'+date+'.root') 
    stop_b = time.time()

    signal_list = process_the_samples(signals, None, None)
    hist_signal = get_histograms(signal_list, variables, None)

    write_hists_to_file(hist_signal, './output_signal_wjet_1l_control_'+date+'.root') 
    stop_s = time.time()

    print "background: ", stop_b - start_b
    print "signal:     ", stop_s - stop_b
    print "total:      ", stop_s - start_b
 
    print 'finished writing'
