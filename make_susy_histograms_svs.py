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
            #hist[sample][tree_name]['N_bjets'] = rt.TH1D('N_bjets_'+sample+'_'+tree_name, 'N S b jets', 20, 0, 20)
            #hist[sample][tree_name]['NS_bjets'] = rt.TH1D('NS_bjets_'+sample+'_'+tree_name, 'N S b jets', 20, 0, 20)
            #hist[sample][tree_name]['NISR_bjets'] = rt.TH1D('NISR_bjets_'+sample+'_'+tree_name, 'N ISR b jets', 20, 0, 20)

            #hist[sample][tree_name]['N_Bs'] = rt.TH1D('N_Bs_'+sample+'_'+tree_name, 'N S Bs', 20, 0, 20)
            #hist[sample][tree_name]['NS_Bs'] = rt.TH1D('NS_Bs_'+sample+'_'+tree_name, 'N S Bs', 20, 0, 20)
            #hist[sample][tree_name]['NISR_Bs'] = rt.TH1D('NISR_Bs_'+sample+'_'+tree_name, 'N ISR Bs', 20, 0, 20)

            hist[sample][tree_name]['Eta_SVs'] = rt.TH1D('Eta_SVs_'+sample+'_'+tree_name, 'N SVs', 64, -3, 3)
            hist[sample][tree_name]['Eta_S_SVs'] = rt.TH1D('Eta_S_SVs_'+sample+'_'+tree_name, 'N S SVs', 64, -3, 3)
            hist[sample][tree_name]['Eta_ISR_SVs'] = rt.TH1D('Eta_ISR_SVs_'+sample+'_'+tree_name, 'N S SVs', 64, -3, 3)
            hist[sample][tree_name]['Sip3DSig'] = rt.TH1D('Sip3DSSig_'+sample+'_'+tree_name, 'N S SVs', 1024, 0, 50)
            
        for ifile, in_file in enumerate(list_of_files_[sample]['files']):
            for tree_name in list_of_files_[sample]['trees']:
                sample_array = get_tree_info_singular(sample, in_file, tree_name, variable_list_, cuts_to_apply_)
                if sample_array is None: continue

                print '\nGetting Histograms for:', sample, tree_name, in_file
                print 'file: ', ifile+1, ' / ', len(list_of_files_[sample]['files'])

                met = np.array(sample_array['MET'])
                risr = np.array(sample_array['RISR'])
                if len(risr) == 0: continue 
                ptisr = np.array(sample_array['PTISR'])

                nlep = np.array(sample_array['Nlep'])
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
                ptisr = np.array([entry[:2] for entry in ptisr])

                nsv_s = np.array([entry[:2] for entry in nsv_s])
                nsv_isr = np.array([entry[:2] for entry in nsv_isr])
                nbjet_s = np.array([entry[:2] for entry in nbjet_s])
                nbjet_isr = np.array([entry[:2] for entry in nbjet_isr])
                isr_index_sv = np.array([entry[:2] if np.shape(entry[:2]) != (0,) else [np.array([],dtype=np.int32),np.array([],dtype=np.int32)] for entry in isr_index_sv])
                s_index_sv = np.array([entry[:2] if np.shape(entry[:2]) != (0,) else [np.array([], dtype=np.int32),np.array([],dtype=np.int32)] for entry in s_index_sv])

                risr = risr[:, 1]
                ptisr = ptisr[:, 1]

                nsv_s = nsv_s[:, 1]
                nsv_isr = nsv_isr[:, 1]
                nbjet_s = nbjet_s[:, 1]
                nbjet_isr = nbjet_isr[:, 1]
                isr_index_sv = isr_index_sv[:, 1]
                s_index_sv = s_index_sv[:, 1]

                weight = 137. * weight
                
                one_lep = nlep == 1 
                two_lep = nlep == 2 
                ################        Medium       #####################
                mini_lep = np.array([mini[lid>=3]*pt[lid>=3] for mini, lid, pt in zip(mini_lep, id_lep, pt_lep)])
                sip3d_lep = np.array([ips[lid>=3] for ips, lid in zip(sip3d_lep, id_lep)])

                #pt_med_lep = np.array([lep[lid>=2] for lep, lid in zip(pt_lep, id_lep)]) 
                #pt_mini_lep = np.array([lep[mini<6] for lep, mini in zip(pt_med_lep, mini_lep)])
                #eta_med_lep = np.array([lep[lid>=2] for lep, lid in zip(eta_lep, id_lep)]) 
                #eta_mini_lep = np.array([lep[mini<6] for lep, mini in zip(eta_med_lep, mini_lep)])
                #phi_med_lep = np.array([lep[lid>=2] for lep, lid in zip(phi_lep, id_lep)]) 
                #phi_mini_lep = np.array([lep[mini<6] for lep, mini in zip(phi_med_lep, mini_lep)])
                #m_med_lep = np.array([lep[lid>=2] for lep, lid in zip(m_lep, id_lep)]) 
                #m_mini_lep = np.array([lep[mini<0.1] for lep, mini in zip(m_med_lep, mini_lep)])


                mini_med_mask = np.array([ True if len(mini[mini<6]) == 1 and one else False for mini, one in zip(mini_lep,one_lep)])
                #mini_med_mask = np.array([ True if len(mini[np.logical_and(mini<6, sip3d<4.)]) == 1 and one else False for mini, sip3d, one in zip(mini_lep,sip3d_lep,one_lep)])
                #mini_med_mask = np.array([ True if len(mini[np.logical_and(mini<6, sip3d<4.)]) == 2 and two else False for mini, sip3d, two in zip(mini_lep,sip3d_lep,two_lep)])
               
 
                eta_sv_s = np.array([eta[ind] for eta, ind in zip(eta_sv, s_index_sv)])
                eta_sv_isr = np.array([eta[ind] for eta, ind in zip(eta_sv, isr_index_sv)])


                risr_0p95 = risr > 0.95
                risr_0p8 = risr > 0.8
                met_200 = met > 200
                ptisr_200 = ptisr > 200
                no_isrb = nbjet_isr < 1
                evt_isrb_selection_mask = np.all([met_200, risr_0p95, ptisr_200, mini_med_mask, no_isrb], axis=0)

                isrb_weight = weight[evt_isrb_selection_mask]

                nbjet = nbjet[evt_isrb_selection_mask]
                nbjet_s = nbjet_s[evt_isrb_selection_mask]
                nbjet_isr = nbjet_isr[evt_isrb_selection_mask]
                nsv_s = nsv_s[evt_isrb_selection_mask]
                nsv = nsv[evt_isrb_selection_mask]
                nsv_isr = nsv_isr[evt_isrb_selection_mask]
                nb = nbjet + nsv
                nb_s = nbjet_s + nsv_s
                nb_isr = nbjet_isr + nsv_isr

                if np.any(evt_isrb_selection_mask):
                    sip3d_lep = np.concatenate(np.array([sip for sip, select in zip(sip3d_lep, evt_isrb_selection_mask) if select]))
                    eta_sv = np.array([eta for eta, select in zip(eta_sv, evt_isrb_selection_mask) if select])
                    eta_sv_s = np.array([eta for eta, select in zip(eta_sv_s, evt_isrb_selection_mask) if select])
                    eta_sv_isr = np.array([eta for eta, select in zip(eta_sv_isr, evt_isrb_selection_mask) if select])
                    isrb_weight_sv = np.array([ [w]*len(sv) for w, sv in zip(isrb_weight, eta_sv)])
                    isrb_weight_s_sv = np.array([ [w]*len(sv) for w, sv in zip(isrb_weight, eta_sv_s)])
                    isrb_weight_isr_sv = np.array([ [w]*len(sv) for w, sv in zip(isrb_weight, eta_sv_isr)])

                    eta_sv = np.concatenate(eta_sv)
                    eta_sv_s = np.concatenate(eta_sv_s)
                    eta_sv_isr = np.concatenate(eta_sv_isr)
                    isrb_weight_sv = np.concatenate(isrb_weight_sv)
                    isrb_weight_s_sv = np.concatenate(isrb_weight_s_sv)
                    isrb_weight_isr_sv = np.concatenate(isrb_weight_isr_sv)

                    rnp.fill_hist(hist[sample][tree_name]['Eta_SVs'], eta_sv, isrb_weight_sv) 
                    rnp.fill_hist(hist[sample][tree_name]['Eta_S_SVs'], eta_sv_s, isrb_weight_s_sv) 
                    rnp.fill_hist(hist[sample][tree_name]['Eta_ISR_SVs'], eta_sv_isr, isrb_weight_isr_sv) 
                    rnp.fill_hist(hist[sample][tree_name]['Sip3DSig'], sip3d_lep, isrb_weight) 
                    #rnp.fill_hist(hist[sample][tree_name]['N_bjets'], nbjet, isrb_weight) 
                    #rnp.fill_hist(hist[sample][tree_name]['NS_bjets'], nbjet_s, isrb_weight) 
                    #rnp.fill_hist(hist[sample][tree_name]['NISR_bjets'], nbjet_isr, isrb_weight) 
                    #rnp.fill_hist(hist[sample][tree_name]['N_SVs'], nsv, isrb_weight) 
                    #rnp.fill_hist(hist[sample][tree_name]['NS_SVs'], nsv_s, isrb_weight) 
                    #rnp.fill_hist(hist[sample][tree_name]['NISR_SVs'], nsv_isr, isrb_weight) 
                    #rnp.fill_hist(hist[sample][tree_name]['N_Bs'], nb, isrb_weight) 
                    #rnp.fill_hist(hist[sample][tree_name]['NS_Bs'], nb_s, isrb_weight) 
                    #rnp.fill_hist(hist[sample][tree_name]['NISR_Bs'], nb_isr, isrb_weight) 

                print 'finished filling'
    return hist


             

if __name__ == "__main__":

    signals = { 
    'SMS-T2bW_dM' : [
                    #'/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_SMS_Stop/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_94X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_102X.root',
], 
    'SMS-T2-4bd_420' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/erichjs/work/Ewkinos/reducer/CMSSW_10_1_4_patch1/src/KUEWKinoAnalysis_dev_v2/output_10Sep19/Fall17_102X_SMS_I/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
],
    'SMS-T2-4bd_490' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/erichjs/work/Ewkinos/reducer/CMSSW_10_1_4_patch1/src/KUEWKinoAnalysis_dev_v2/output_10Sep19/Fall17_102X_SMS_I/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
],
#    'SMS-TChiWH' : [ 
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/SMS-TChiWH_WToLNu_HToBB_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_102X',
#],
#    'SMS-T2tt_dM' : [
#                    #'/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/Fall17_102X_SMS_Stop/SMS-T2tt_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/SMS-T2tt_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_102X.root',
#], 
              }
    backgrounds = {
    'TTJets_2017' : [
#                     '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X'
#                     '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X'
#                     '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/Fall17_102X_bkgextra/TTTT_TuneCP5_PSweights_13TeV-amcatnlo-pythia8_Fall17_102X',
                     '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                     '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                     '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
],
    'WJets_2017' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/NEW_29_11_19/NoHadd/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
],
#    'ZJets_2017' : [
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-100To200_13TeV-madgraph_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-200To400_13TeV-madgraph_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-400To600_13TeV-madgraph_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-600To800_13TeV-madgraph_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-800To1200_13TeV-madgraph_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph_Fall17_94X',
#],
#    'DY_M50_2017' : [
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_94X',
#],
#    'WW_2017' : [
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WWTo4Q_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WWG_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_94X',
#],
#    'ZZ_2017' : [
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZZTo2L2Nu_13TeV_powheg_pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZZTo2Q2Nu_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZZTo4L_13TeV_powheg_pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_94X',
#],
#    'WZ_2017' : [
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZG_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZZ_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_94X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_94X_bkg/WZ_TuneCP5_13TeV-pythia8_Fall17_94X',
#],
                  }
    variables = ['MET', 'RISR', 'PTISR', 'PT_lep', 'MiniIso_lep', 'ID_lep', 'Nbjet', 'Nbjet_S', 'Nbjet_ISR', 'NSV_S', 'NSV_ISR', 'weight', 'NSV', 'Nlep', 'Eta_SV', 'index_SV_S','index_SV_ISR', 'SIP3D_lep']

    start_b = time.time()    
    background_list = process_the_samples(backgrounds, None, None)
    hist_background = get_histograms(background_list, variables, None)

    write_hists_to_file(hist_background, './output_background_risr_0p8_mixed_eta_'+date+'.root') 
    stop_b = time.time()

    signal_list = process_the_samples(signals, None, None)
    hist_signal = get_histograms(signal_list, variables, None)

    write_hists_to_file(hist_signal, './output_signal_risr_0p8_mixed_eta_'+date+'.root') 
    stop_s = time.time()

    print "background: ", stop_b - start_b
    print "signal:     ", stop_s - stop_b
    print "total:      ", stop_s - start_b
 
    print 'finished writing'
