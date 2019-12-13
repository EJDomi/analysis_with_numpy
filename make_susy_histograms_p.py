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
            hist[sample][tree_name]['PS'] = rt.TH1D('PS_'+sample+'_'+tree_name, 'P S system', 1024, 0, 500)
            hist[sample][tree_name]['PV'] = rt.TH1D('PS_'+sample+'_'+tree_name, 'P S system', 1024, 0, 500)
            hist[sample][tree_name]['P_lep'] = rt.TH1D('P_lep_'+sample+'_'+tree_name, 'P leps', 1024, 0, 500)
            hist[sample][tree_name]['Eta_SVs'] = rt.TH1D('Eta_SVs_'+sample+'_'+tree_name, 'N SVs', 64, -3, 3)
            hist[sample][tree_name]['Eta_SVs_prerisr'] = rt.TH1D('Eta_SVs_prerisr_'+sample+'_'+tree_name, 'N SVs', 64, -3, 3)
            hist[sample][tree_name]['Eta_SVs_preisrb'] = rt.TH1D('Eta_SVs_preisrb'+sample+'_'+tree_name, 'N SVs', 64, -3, 3)
            
        for ifile, in_file in enumerate(list_of_files_[sample]['files']):
            for tree_name in list_of_files_[sample]['trees']:
                sample_array = get_tree_info_singular(sample, in_file, tree_name, variable_list_, cuts_to_apply_)
                if sample_array is None: continue

                print '\nGetting Histograms for:', sample, tree_name, in_file
                print 'file: ', ifile+1, ' / ', len(list_of_files_[sample]['files'])

                met = np.array(sample_array['MET'])
                ps = np.array(sample_array['PS'])
                pv = np.array(sample_array['PV'])
                pla = np.array(sample_array['PLa'])
                plb = np.array(sample_array['PLb'])
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

                weight = np.array(sample_array['weight'])

                id_lep = np.array(sample_array['ID_lep'])
                pt_lep = np.array(sample_array['PT_lep'])
                mini_lep = np.array(sample_array['MiniIso_lep'])
                sip3d_lep = np.array(sample_array['SIP3D_lep'])
                
                risr = np.array([entry[:2] for entry in risr])
                ps = np.array([entry[:2] for entry in ps])
                pv = np.array([entry[:2] for entry in pv])
                pla = np.array([entry[:2] for entry in pla])
                plb = np.array([entry[:2] for entry in plb])
                ptisr = np.array([entry[:2] for entry in ptisr])

                nsv_s = np.array([entry[:2] for entry in nsv_s])
                nsv_isr = np.array([entry[:2] for entry in nsv_isr])
                nbjet_s = np.array([entry[:2] for entry in nbjet_s])
                nbjet_isr = np.array([entry[:2] for entry in nbjet_isr])

                risr = risr[:, 1]
                ps = ps[:, 1]
                pv = pv[:, 1]
                pla = pla[:, 1]
                plb = plb[:, 1]
                ptisr = ptisr[:, 1]

                nsv_s = nsv_s[:, 1]
                nsv_isr = nsv_isr[:, 1]
                nbjet_s = nbjet_s[:, 1]
                nbjet_isr = nbjet_isr[:, 1]

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
               
 
                risr_0p95 = risr > 0.95
                risr_0p8 = risr > 0.8
                met_200 = met > 200
                ptisr_200 = ptisr > 200
                no_isrb = nbjet_isr < 1
                evt_isrb_selection_mask = np.all([met_200, risr_0p95, ptisr_200, mini_med_mask, no_isrb], axis=0)
                evt_preisrb_selection_mask = np.all([met_200, risr_0p95, ptisr_200, mini_med_mask], axis=0)
                evt_prerisr_selection_mask = np.all([met_200, ptisr_200, mini_med_mask], axis=0)

                isrb_weight = weight[evt_isrb_selection_mask]
                preisrb_weight = weight[evt_preisrb_selection_mask]
                prerisr_weight = weight[evt_prerisr_selection_mask]

                if np.any(evt_preisrb_selection_mask):
                    eta_preisrb_sv = np.array([eta for eta, select in zip(eta_sv, evt_preisrb_selection_mask) if select])
                    preisrb_weight_sv = np.array([ [w]*len(sv) for w, sv in zip(preisrb_weight, eta_preisrb_sv)])

                    eta_preisrb_sv = np.concatenate(eta_preisrb_sv)
                    preisrb_weight_sv = np.concatenate(preisrb_weight_sv)

                    rnp.fill_hist(hist[sample][tree_name]['Eta_SVs_preisrb'], eta_preisrb_sv, preisrb_weight_sv) 

                if np.any(evt_prerisr_selection_mask):
                    eta_prerisr_sv = np.array([eta for eta, select in zip(eta_sv, evt_prerisr_selection_mask) if select])
                    prerisr_weight_sv = np.array([ [w]*len(sv) for w, sv in zip(prerisr_weight, eta_prerisr_sv)])

                    eta_prerisr_sv = np.concatenate(eta_prerisr_sv)
                    prerisr_weight_sv = np.concatenate(prerisr_weight_sv)

                    rnp.fill_hist(hist[sample][tree_name]['Eta_SVs_prerisr'], eta_prerisr_sv, prerisr_weight_sv)
 
                if np.any(evt_isrb_selection_mask):
                    eta_all_sv = np.array([eta for eta, select in zip(eta_sv, evt_isrb_selection_mask) if select])
                    isrb_weight_sv = np.array([ [w]*len(sv) for w, sv in zip(isrb_weight, eta_all_sv)])

                    ps = ps[evt_isrb_selection_mask]
                    pv = pv[evt_isrb_selection_mask]
                    pla = pla[evt_isrb_selection_mask]
                    plb = plb[evt_isrb_selection_mask]

                    pl = pla + plb

                    eta_all_sv = np.concatenate(eta_all_sv)
                    isrb_weight_sv = np.concatenate(isrb_weight_sv)
                    print len(pla), len(plb), len(pl)
                    rnp.fill_hist(hist[sample][tree_name]['Eta_SVs'], eta_all_sv, isrb_weight_sv) 
                    rnp.fill_hist(hist[sample][tree_name]['PS'], ps, isrb_weight) 
                    rnp.fill_hist(hist[sample][tree_name]['PV'], pv, isrb_weight) 
                    rnp.fill_hist(hist[sample][tree_name]['P_lep'], pl, isrb_weight) 

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
    variables = ['MET', 'RISR', 'PTISR', 'PT_lep', 'MiniIso_lep', 'ID_lep', 'Nbjet', 'Nbjet_S', 'Nbjet_ISR', 'NSV_S', 'NSV_ISR', 'weight', 'NSV', 'Nlep', 'Eta_SV', 'index_SV_S','index_SV_ISR', 'SIP3D_lep', 'PS', 'PV', 'PLa', 'PLb']

    start_b = time.time()    
    background_list = process_the_samples(backgrounds, None, None)
    hist_background = get_histograms(background_list, variables, None)

    write_hists_to_file(hist_background, './output_background_risr_0p8_mixed_p_'+date+'.root') 
    stop_b = time.time()

    signal_list = process_the_samples(signals, None, None)
    hist_signal = get_histograms(signal_list, variables, None)

    write_hists_to_file(hist_signal, './output_signal_risr_0p8_mixed_p_'+date+'.root') 
    stop_s = time.time()

    print "background: ", stop_b - start_b
    print "signal:     ", stop_s - stop_b
    print "total:      ", stop_s - start_b
 
    print 'finished writing'
