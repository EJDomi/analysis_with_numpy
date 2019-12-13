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
            hist[sample][tree_name]['PT_b'] = rt.TH1D('PT_b_'+sample+'_'+tree_name, 'p_{T} gen bs', 512, 0, 1000)
            hist[sample][tree_name]['PT_W'] = rt.TH1D('PT_W_'+sample+'_'+tree_name, 'p_{T} gen Ws', 512, 0, 1000)
            hist[sample][tree_name]['PT_chi'] = rt.TH1D('PT_chi_'+sample+'_'+tree_name, 'p_{T} gen chis', 512, 0, 1000)
            
        for ifile, in_file in enumerate(list_of_files_[sample]['files']):
            for tree_name in list_of_files_[sample]['trees']:
                sample_array = get_tree_info_singular(sample, in_file, tree_name, variable_list_, cuts_to_apply_)
                if sample_array is None: continue

                print '\nGetting Histograms for:', sample, tree_name, in_file
                print 'file: ', ifile+1, ' / ', len(list_of_files_[sample]['files'])

                gen_pdgid = np.array(sample_array['recoGenParticles_prunedGenParticles__PAT.obj.m_state.pdgId_'])
                gen_pt = np.array(sample_array['recoGenParticles_prunedGenParticles__PAT.obj.m_state.p4Polar_.fCoordinates.fPt'])
                gen_eta = np.array(sample_array['recoGenParticles_prunedGenParticles__PAT.obj.m_state.p4Polar_.fCoordinates.fEta'])
                gen_mass = np.array(sample_array['recoGenParticles_prunedGenParticles__PAT.obj.m_state.p4Polar_.fCoordinates.fM'])

                is_420_susy = np.full(len(gen_pt), False)
                is_460_susy = np.full(len(gen_pt), False)
                is_490_susy = np.full(len(gen_pt), False)

                if 'SMS_T2bW' in sample:
                    stop_mass = np.array([mass[np.logical_or(pdgid == 1000006, pdgid == 2000006)] for mass, pdgid in zip(gen_mass, gen_pdgid)])
                    child_mass = np.array([mass[pdgid == 1000022] for mass, pdgid in zip(gen_mass, gen_pdgid)])
                    is_420_susy = np.array([True if np.any(stop == 500.) and np.any(child == 420.) else False for stop, child in zip(stop_mass, child_mass)])
                    is_460_susy = np.array([True if np.any(stop == 500.) and np.any(child == 460.) else False for stop, child in zip(stop_mass, child_mass)])
                    is_490_susy = np.array([True if np.any(stop == 500.) and np.any(child == 490.) else False for stop, child in zip(stop_mass, child_mass)])
                
 

                if np.any(is_420_susy) and 'T2bW_420' in sample:
                    b_pt = np.array([pt[np.logical_and(np.abs(pdgid) == 5, np.abs(eta)<2.4)] for pt, eta, pdgid, susy in zip(gen_pt, gen_eta, gen_pdgid, is_420_susy) if susy])
                    chi_pt = np.array([pt[np.logical_and(np.abs(pdgid) == 5, np.abs(eta)<2.4)] for pt, eta, pdgid, susy in zip(gen_pt, gen_eta, gen_pdgid, is_420_susy) if susy])
                    w_pt = np.array([pt[np.logical_and(np.abs(pdgid) == 24, np.abs(eta)<2.4)] for pt, eta, pdgid, susy in zip(gen_pt, gen_eta, gen_pdgid, is_420_susy) if susy])
                    rnp.fill_hist(hist[sample][tree_name]['PT_b'], np.concatenate(b_pt))
                    rnp.fill_hist(hist[sample][tree_name]['PT_W'], np.concatenate(w_pt))
                    rnp.fill_hist(hist[sample][tree_name]['PT_chi'], np.concatenate(chi_pt))
                elif np.any(is_460_susy) and 'T2bW_460' in sample:
                    b_pt = np.array([pt[np.logical_and(np.abs(pdgid) == 5, np.abs(eta)<2.4)] for pt, eta, pdgid, susy in zip(gen_pt, gen_eta, gen_pdgid, is_460_susy) if susy])
                    chi_pt = np.array([pt[np.logical_and(np.abs(pdgid) == 5, np.abs(eta)<2.4)] for pt, eta, pdgid, susy in zip(gen_pt, gen_eta, gen_pdgid, is_460_susy) if susy])
                    w_pt = np.array([pt[np.logical_and(np.abs(pdgid) == 24, np.abs(eta)<2.4)] for pt, eta, pdgid, susy in zip(gen_pt, gen_eta, gen_pdgid, is_460_susy) if susy])
                    rnp.fill_hist(hist[sample][tree_name]['PT_b'], np.concatenate(b_pt))
                    rnp.fill_hist(hist[sample][tree_name]['PT_W'], np.concatenate(w_pt))
                    rnp.fill_hist(hist[sample][tree_name]['PT_chi'], np.concatenate(chi_pt))
                elif np.any(is_490_susy) and 'T2bW_490' in sample:
                    b_pt = np.array([pt[np.logical_and(np.abs(pdgid) == 5, np.abs(eta)<2.4)] for pt, eta, pdgid, susy in zip(gen_pt, gen_eta, gen_pdgid, is_490_susy) if susy])
                    chi_pt = np.array([pt[np.logical_and(np.abs(pdgid) == 5, np.abs(eta)<2.4)] for pt, eta, pdgid, susy in zip(gen_pt, gen_eta, gen_pdgid, is_490_susy) if susy])
                    w_pt = np.array([pt[np.logical_and(np.abs(pdgid) == 24, np.abs(eta)<2.4)] for pt, eta, pdgid, susy in zip(gen_pt, gen_eta, gen_pdgid, is_490_susy) if susy])
                    rnp.fill_hist(hist[sample][tree_name]['PT_b'], np.concatenate(b_pt))
                    rnp.fill_hist(hist[sample][tree_name]['PT_W'], np.concatenate(w_pt))
                    rnp.fill_hist(hist[sample][tree_name]['PT_chi'], np.concatenate(chi_pt))
                else:
                    b_pt = np.array([pt[np.logical_and(np.abs(pdgid) == 5, np.abs(eta)<2.4)] for pt, eta, pdgid in zip(gen_pt, gen_eta, gen_pdgid)])
                    if 'SMS' in sample:
                        chi_pt = np.array([pt[np.logical_and(pdgid == 1000022, np.abs(eta)<2.4)] for pt, eta, pdgid in zip(gen_pt, gen_eta, gen_pdgid)])
                        rnp.fill_hist(hist[sample][tree_name]['PT_chi'], np.concatenate(chi_pt))
                    w_pt = np.array([pt[np.logical_and(np.abs(pdgid) == 24, np.abs(eta)<2.4)] for pt, eta, pdgid in zip(gen_pt, gen_eta, gen_pdgid)])
                    rnp.fill_hist(hist[sample][tree_name]['PT_b'], np.concatenate(b_pt))
                    rnp.fill_hist(hist[sample][tree_name]['PT_W'], np.concatenate(w_pt))

                print 'finished filling'
    return hist


             

if __name__ == "__main__":

    signals = { 
    'SMS_T2bW_420' : [
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/FE75BCEF-6B47-E911-9021-0CC47A7C347E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/FAF80FCB-4F47-E911-A04D-848F69FD0B51.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/FAD5BE4C-8446-E911-B043-0CC47A7C3408.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/F6E18E0C-6A46-E911-A7B9-E0071B7B2350.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/F4735C5C-8A46-E911-8988-6C3BE5B59058.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/F2342584-6B47-E911-955A-509A4C74D1B3.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/EE6CC173-F746-E911-BE96-44A84225C893.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E80D5CA7-7A46-E911-961E-0CC47AFB7DCC.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E6F584CF-6547-E911-8A39-AC1F6B0DE33A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E6681DF6-B146-E911-808F-001F29089F7E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E28D2518-7146-E911-A396-FA163E7DC8D8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E0FA8C44-5847-E911-BF78-90B11CBD0004.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E002A7E8-3747-E911-AD25-0CC47A7DFEDA.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/DE3E9654-4D47-E911-B4C0-0242AC1C0501.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/DC43D654-6F46-E911-9F9F-FA163E0734FD.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/DA00EC89-7346-E911-B5D5-40F2E9C6AE21.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/D691F22F-7646-E911-B5FB-FA163E667F65.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/D6255F44-0B47-E911-A889-20CF305B064B.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/D0C3322B-7846-E911-B144-B8CA3A709648.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/CEAFEF46-7E46-E911-BE0A-A0369FC5E0A4.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/CAFA126D-7546-E911-B5D9-02163E01A01A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/CA9EF27E-6C46-E911-9C8E-24BE05CEECD1.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/C8841A7A-8D46-E911-B349-C45444922D6C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/C626EB52-6B47-E911-B030-AC1F6B1AF05A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/BE49136F-6C46-E911-8A8A-FA163E1B2960.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/B6C1DF0C-9746-E911-A4A0-001E6779258C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/B4C5D7B2-6046-E911-A648-FA163EF36544.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/B48CD991-7C46-E911-9B90-008CFAC93CB0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/AC9B42A2-5A47-E911-9FBE-1CB72C1B6C46.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/AC4A1F1F-7B46-E911-BB88-0025905D1E08.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/A8C4C047-5847-E911-84F4-0025904C7F82.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/A6826D99-FD46-E911-A9A5-A0369FD0B35A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/A4218663-6F46-E911-AB84-FA163ECB6712.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9E0DADBE-6646-E911-85E9-FA163E29485B.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9AC812A4-7846-E911-9DD0-FA163EC9203A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9A975843-6A46-E911-BABE-FA163EAF5860.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/964F51F8-7B46-E911-8928-002590E3A0FA.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/94F04AF1-4847-E911-B949-008CFA197448.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/944B84D5-7B46-E911-9994-FA163E85491A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/92A4550D-4347-E911-86E3-0CC47AFF02F8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9096F370-6646-E911-A8BE-FA163EA40BD4.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/90380585-6B47-E911-B62C-AC1F6B8DBEC2.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/88368126-8146-E911-9DB1-0CC47A4C8F30.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/84EDC0C3-8746-E911-BC9C-A4BF0112DFA0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/7866B56F-7546-E911-A5E2-0242AC1C0503.root',
],
    'SMS_T2bW_460' : [
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/FE75BCEF-6B47-E911-9021-0CC47A7C347E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/FAF80FCB-4F47-E911-A04D-848F69FD0B51.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/FAD5BE4C-8446-E911-B043-0CC47A7C3408.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/F6E18E0C-6A46-E911-A7B9-E0071B7B2350.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/F4735C5C-8A46-E911-8988-6C3BE5B59058.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/F2342584-6B47-E911-955A-509A4C74D1B3.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/EE6CC173-F746-E911-BE96-44A84225C893.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E80D5CA7-7A46-E911-961E-0CC47AFB7DCC.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E6F584CF-6547-E911-8A39-AC1F6B0DE33A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E6681DF6-B146-E911-808F-001F29089F7E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E28D2518-7146-E911-A396-FA163E7DC8D8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E0FA8C44-5847-E911-BF78-90B11CBD0004.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E002A7E8-3747-E911-AD25-0CC47A7DFEDA.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/DE3E9654-4D47-E911-B4C0-0242AC1C0501.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/DC43D654-6F46-E911-9F9F-FA163E0734FD.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/DA00EC89-7346-E911-B5D5-40F2E9C6AE21.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/D691F22F-7646-E911-B5FB-FA163E667F65.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/D6255F44-0B47-E911-A889-20CF305B064B.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/D0C3322B-7846-E911-B144-B8CA3A709648.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/CEAFEF46-7E46-E911-BE0A-A0369FC5E0A4.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/CAFA126D-7546-E911-B5D9-02163E01A01A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/CA9EF27E-6C46-E911-9C8E-24BE05CEECD1.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/C8841A7A-8D46-E911-B349-C45444922D6C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/C626EB52-6B47-E911-B030-AC1F6B1AF05A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/BE49136F-6C46-E911-8A8A-FA163E1B2960.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/B6C1DF0C-9746-E911-A4A0-001E6779258C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/B4C5D7B2-6046-E911-A648-FA163EF36544.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/B48CD991-7C46-E911-9B90-008CFAC93CB0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/AC9B42A2-5A47-E911-9FBE-1CB72C1B6C46.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/AC4A1F1F-7B46-E911-BB88-0025905D1E08.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/A8C4C047-5847-E911-84F4-0025904C7F82.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/A6826D99-FD46-E911-A9A5-A0369FD0B35A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/A4218663-6F46-E911-AB84-FA163ECB6712.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9E0DADBE-6646-E911-85E9-FA163E29485B.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9AC812A4-7846-E911-9DD0-FA163EC9203A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9A975843-6A46-E911-BABE-FA163EAF5860.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/964F51F8-7B46-E911-8928-002590E3A0FA.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/94F04AF1-4847-E911-B949-008CFA197448.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/944B84D5-7B46-E911-9994-FA163E85491A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/92A4550D-4347-E911-86E3-0CC47AFF02F8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9096F370-6646-E911-A8BE-FA163EA40BD4.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/90380585-6B47-E911-B62C-AC1F6B8DBEC2.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/88368126-8146-E911-9DB1-0CC47A4C8F30.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/84EDC0C3-8746-E911-BC9C-A4BF0112DFA0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/7866B56F-7546-E911-A5E2-0242AC1C0503.root',
],
    'SMS_T2bW_490' : [
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/94F04AF1-4847-E911-B949-008CFA197448.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/944B84D5-7B46-E911-9994-FA163E85491A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/92A4550D-4347-E911-86E3-0CC47AFF02F8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9096F370-6646-E911-A8BE-FA163EA40BD4.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/90380585-6B47-E911-B62C-AC1F6B8DBEC2.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/88368126-8146-E911-9DB1-0CC47A4C8F30.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/84EDC0C3-8746-E911-BC9C-A4BF0112DFA0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/7866B56F-7546-E911-A5E2-0242AC1C0503.root',
],
    'SMS_T2bW_460' : [
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/FE75BCEF-6B47-E911-9021-0CC47A7C347E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/FAF80FCB-4F47-E911-A04D-848F69FD0B51.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/FAD5BE4C-8446-E911-B043-0CC47A7C3408.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/F6E18E0C-6A46-E911-A7B9-E0071B7B2350.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/F4735C5C-8A46-E911-8988-6C3BE5B59058.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/F2342584-6B47-E911-955A-509A4C74D1B3.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/EE6CC173-F746-E911-BE96-44A84225C893.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E80D5CA7-7A46-E911-961E-0CC47AFB7DCC.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E6F584CF-6547-E911-8A39-AC1F6B0DE33A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E6681DF6-B146-E911-808F-001F29089F7E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E28D2518-7146-E911-A396-FA163E7DC8D8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E0FA8C44-5847-E911-BF78-90B11CBD0004.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/E002A7E8-3747-E911-AD25-0CC47A7DFEDA.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/DE3E9654-4D47-E911-B4C0-0242AC1C0501.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/DC43D654-6F46-E911-9F9F-FA163E0734FD.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/DA00EC89-7346-E911-B5D5-40F2E9C6AE21.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/D691F22F-7646-E911-B5FB-FA163E667F65.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/D6255F44-0B47-E911-A889-20CF305B064B.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/D0C3322B-7846-E911-B144-B8CA3A709648.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/CEAFEF46-7E46-E911-BE0A-A0369FC5E0A4.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/CAFA126D-7546-E911-B5D9-02163E01A01A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/CA9EF27E-6C46-E911-9C8E-24BE05CEECD1.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/C8841A7A-8D46-E911-B349-C45444922D6C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/C626EB52-6B47-E911-B030-AC1F6B1AF05A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/BE49136F-6C46-E911-8A8A-FA163E1B2960.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/B6C1DF0C-9746-E911-A4A0-001E6779258C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/B4C5D7B2-6046-E911-A648-FA163EF36544.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/B48CD991-7C46-E911-9B90-008CFAC93CB0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/AC9B42A2-5A47-E911-9FBE-1CB72C1B6C46.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/AC4A1F1F-7B46-E911-BB88-0025905D1E08.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/A8C4C047-5847-E911-84F4-0025904C7F82.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/A6826D99-FD46-E911-A9A5-A0369FD0B35A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/A4218663-6F46-E911-AB84-FA163ECB6712.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9E0DADBE-6646-E911-85E9-FA163E29485B.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9AC812A4-7846-E911-9DD0-FA163EC9203A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9A975843-6A46-E911-BABE-FA163EAF5860.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/964F51F8-7B46-E911-8928-002590E3A0FA.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/94F04AF1-4847-E911-B949-008CFA197448.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/944B84D5-7B46-E911-9994-FA163E85491A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/92A4550D-4347-E911-86E3-0CC47AFF02F8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/9096F370-6646-E911-A8BE-FA163EA40BD4.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/90380585-6B47-E911-B62C-AC1F6B8DBEC2.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/88368126-8146-E911-9DB1-0CC47A4C8F30.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/84EDC0C3-8746-E911-BC9C-A4BF0112DFA0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUFall17Fast_94X_mc2017_realistic_v15-v1/60000/7866B56F-7546-E911-A5E2-0242AC1C0503.root',
],
    'SMS_T2_4bd_420' : [
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/E28F8FEE-2EA6-E811-86DA-A4BF0112BCBA.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/DA29F869-CBA7-E811-885C-FA163EE8F7A6.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/94DFDF9F-CAA7-E811-8DD6-002590DE6E2E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/8007411A-CCA7-E811-B4AB-AC1F6B1AF05A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/6C97ED03-8BA4-E811-9163-A4BF0112DC34.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/247C6C51-64A8-E811-B235-008CFA580778.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/1603F4EA-65A8-E811-ABE0-14187741278B.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/F6CA5F8D-E79D-E811-ABF7-A4BF01125730.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/E4EECE59-DE9F-E811-823D-3417EBE70729.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/E257B132-C39D-E811-A049-00266CFFBC74.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/DA2FA442-789D-E811-861D-0025905C5484.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/D08701B7-A99C-E811-B2D8-008CFAFBEEE6.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/C83896BE-DE9F-E811-84BA-008CFA1974DC.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/C2D23CCF-DE9F-E811-8D63-782BCB3BCE9F.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/C0196058-DE9F-E811-B1E7-00259073E522.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B6F7F686-B19C-E811-9E81-001E67E6F8AF.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B46E01AF-DE9F-E811-964B-7CD30AB04506.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/A693BDC2-4E9E-E811-8BA5-0CC47AD98D6C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/9C740EB7-2C9F-E811-A409-0025B3268576.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/90B19EA2-DE9F-E811-9CCD-141877410E71.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/8AF2EDC7-AD9E-E811-A661-3417EBE64BAF.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/8057BDE5-659E-E811-9407-F01FAFE5CC96.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/7C7F5456-4F9D-E811-B3C5-008CFAF29284.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6EAFBE69-E79C-E811-950E-009C029C120A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6AED948E-709D-E811-A145-0CC47AD98F6A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6469BB6D-DE9F-E811-A940-3417EBE669D4.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/56F5C050-749E-E811-B848-1418774121A1.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/4CC36B76-879C-E811-AB72-44A842B2990B.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/4AF4D680-B89E-E811-8D82-A4BF0112BE32.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/48102D32-E99C-E811-B7F7-0026B94DBE17.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/40CC80B9-DE9F-E811-9B18-002590E7E01A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/3ADA79C2-DE9F-E811-8728-EC0D9A0B30D0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/3A12C3AF-DE9F-E811-AE18-B4969109FA60.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/388FD26B-DE9F-E811-902A-0025905C5502.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/3485613B-F49D-E811-B9F7-E0071B693B41.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/30AEA16C-DE9F-E811-9922-AC1F6B1E2F64.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/22CAA4B9-F69C-E811-BD2B-00266CFEFC5C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/207CDF54-479E-E811-8AC1-0025905C975C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/1EC2B083-DE9F-E811-B72A-00266CF89604.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/1E6A8988-0A9E-E811-BE82-D48564597C70.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/149C7670-DE9F-E811-9185-089E01DC6B5C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/145F778B-939D-E811-AC5E-0CC47A0109A6.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/1290A448-659E-E811-9745-0CC47A4D7666.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/12339972-DE9F-E811-9B20-0CC47A7C3612.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/10441A89-839F-E811-BE01-509A4C84EABE.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/08620262-DE9F-E811-9430-90E2BA034674.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/00A888AC-DE9F-E811-9BF4-B8CA3A708F98.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/004C1F72-DE9F-E811-83D3-44A842CFC9B2.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/100000/C83CEA2F-2FB8-E811-A22E-001E6779250C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/100000/701BAC0F-25B8-E811-B469-A4BF01125D8E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/100000/58AE27D6-02B9-E811-A57E-0242AC1C0502.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/100000/1A7C2F95-7AB8-E811-B1A6-A4BF01158960.root',
],
    'SMS_T2_4bd_490' : [
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/F4564557-27A2-E811-A656-0CC47A7DFFA6.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/D88495FC-6AAC-E811-AC31-24BE05C3EC61.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/D26D09D3-6AAC-E811-AC0E-0CC47A5FC67D.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/CAFD8051-A1A7-E811-96E4-FA163E857ACE.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/CAE6DBE9-6AAC-E811-923B-008CFA0514E0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/C61EF366-5A95-E811-9928-FA163E4A464D.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/C2166FD4-6AAC-E811-952A-A0369F5BD91C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/B692D644-75A0-E811-8973-C454449229EB.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/AE7F2D34-B891-E811-880A-0CC47A7E6B14.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/A44478F1-9492-E811-93A2-0017A4770450.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/A08E8054-7795-E811-806D-C81F66B7EBF4.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/9054BDCA-6AAC-E811-8D66-0CC47A7FC702.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/8E9510CF-6AAC-E811-92FD-0025905A60D6.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/8047F748-C48F-E811-83C7-00266CF3E174.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/6C322185-6AAC-E811-957F-001E67792514.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/68E79DD9-6AAC-E811-BD83-44A842CFD626.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/5E22428F-A9A1-E811-8488-008CFA111200.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/52E95CA7-A391-E811-8AA7-0025904C7DF6.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/4E8F0C0E-6CAC-E811-AA59-509A4C7489AF.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/4A4654D3-6AAC-E811-9EB6-002590E7D7EA.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/486A271D-DCA6-E811-813E-001E677923AE.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/38B321C6-F795-E811-8B7F-549F351F915E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/3416468C-6AAC-E811-A071-1866DAEEB0C0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/22179D90-FA90-E811-9A5E-FA163E4AB21B.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/18BF04CD-6AAC-E811-A180-44A84225C7BB.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/105E9D1A-6BAC-E811-9440-FA163EFF440B.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/0A6EE1EF-6AAC-E811-85A9-0CC47AFC3D34.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/086C3B7B-8FA6-E811-A932-0025904CDDEE.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/FAFAF3F5-B0C6-E811-AB3E-A0369F301924.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/FA0E5A71-E8C6-E811-A09E-001E67E6F7C4.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/F28233BD-5EC9-E811-A095-0CC47A78A496.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/DC677FA5-03C7-E811-8253-001E677926FE.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/D0686FA8-03C6-E811-8150-AC1F6B8DD244.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/8E5CE2F7-5EC9-E811-BA80-003048F5ADF0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/8CF08937-A1C7-E811-BADD-0CC47AA9943A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/7A5C2FE7-5EC9-E811-8939-D4AE526A0B47.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/7866BF3F-B6C7-E811-B5CB-A0369F310120.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/7471CC90-EDC6-E811-901F-A4BF0112BDA8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/6A719AE0-5EC9-E811-81B2-24BE05BD4F61.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/5062D3D6-B0C6-E811-A136-00266CF85940.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/44731C2A-D9C5-E811-90D6-3417EBE47FE5.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/4449BF08-A2C8-E811-8594-001E677923A0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/444720EE-58C8-E811-8050-1866DA8797A4.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/361ECE07-B1C6-E811-BEA3-A4BF01025C02.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/245D62E6-5EC9-E811-9A38-001E67504255.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/1E9664D0-B0C6-E811-9880-A0369FD20740.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/08FBB9BC-5EC9-E811-9BBF-A4BF011255F8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/08E83F55-D8C6-E811-95A3-001E67E6F922.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/100000/7E6F39AE-3DC6-E811-9BFD-0025905B85DE.root',
],
              }
    backgrounds = {
    'TTJets_2017' : [
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/FEC6063F-BE88-E811-9D5F-1866DAEB5C78.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/F6E4BBB7-7088-E811-97CB-F01FAFE37F53.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/F6C0DB22-0A88-E811-BED6-0CC47A7C34C4.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/F6B9538E-6B86-E811-B026-0242AC1C0500.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/F483D01D-5889-E811-915F-0CC47A7AB7A0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/F2D0F584-EA88-E811-9B7A-0025905A60A0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/ECA23218-4687-E811-89F0-44A842BE8F98.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/EC4DC386-C292-E811-8385-00266CFFBE5C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/E623A85A-108A-E811-B7D0-D4AE52E94750.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/E620D45B-0B92-E811-B90A-901B0E5427A6.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/E4A4D33B-C292-E811-BE13-44A842CFC9E6.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/E45254F2-A088-E811-A057-0CC47A7C35C8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/E408683D-FB87-E811-8F8D-0025907D2446.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/E0EAD899-C292-E811-9884-14DDA9243247.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/E00E6130-5E87-E811-A787-002590D9D8BA.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/DE0F9C3D-2187-E811-8860-0CC47AD98C5E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/DC59D8B6-AF87-E811-B916-4C79BA181225.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/D6FCD4F3-D58A-E811-ABB3-0242AC1C0501.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/D6BCC468-AF87-E811-A723-001CC445D6D2.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/D45302C2-FB87-E811-8C43-F01FAFE15857.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/D2A4CF1F-C292-E811-A4B6-008CFA111220.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/D0CF3CA2-2387-E811-8530-0242AC130002.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/CEEFB3F9-5789-E811-A44C-A4BF0112BC76.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/C45B0562-E986-E811-86F2-0CC47AD98BEA.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/C2914E69-5889-E811-9562-0CC47A13CCFC.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/C075AFDE-FB87-E811-9D5C-000101000965.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/BEFF8D30-2187-E811-B1DA-1866DAEA8394.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/BCB9A1A7-3387-E811-9136-00010100097A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/BC2CF28A-E286-E811-A7C3-0CC47A7FC6F8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/BAC2E7D2-FB87-E811-B0E7-000101000972.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/B625B999-BB8F-E811-AF42-06FDAA0000B2.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/B4AC32A0-9489-E811-9DCE-24BE05C4D851.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/B05FADA3-C292-E811-9E31-842B2B6FE5EF.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/B03D1CBC-4186-E811-A9A5-24BE05C46B11.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/ACF29CDD-FB87-E811-BAE4-000101000938.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/ACEE07F5-1789-E811-95B1-0CC47A4C8E56.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/AA13E100-9590-E811-8E41-A0369FE2C05E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/A8B9962E-A288-E811-8DC8-1866DA87B230.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/A83124B8-C292-E811-905D-0CC47AD9901E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/A6BF2270-4186-E811-8208-44A842B2D631.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/A6B5F128-5987-E811-83E8-0CC47A57CD56.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/A6A95B84-5987-E811-A53F-0CC47AD98CF6.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/A4EA13A0-D386-E811-A152-0025905C53B2.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/A069B6F6-A586-E811-97D1-003048947BBB.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/9E883389-5889-E811-95A7-00010100096A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/9AD40E10-C392-E811-B92A-001E6739A959.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/98AA966E-AF87-E811-AAA4-0242AC130002.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/96B2C1D7-FB87-E811-AEA0-00010100096F.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/967ED8AB-6387-E811-AD53-008CFA1974A4.root',
],
    'WJets_2017' : [
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/F6612338-BC59-E911-9FFD-0CC47AD98D6E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/F8D2C919-165A-E911-91B0-0CC47A4C8E66.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/FEF1EBDA-9059-E911-B71A-008CFAE452D0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/FE956420-925A-E911-8872-0025905A48FC.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/FE90EFF0-7C5A-E911-BA9D-AC1F6BAC7D1A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/FCB43385-275B-E911-86FE-A0369FE2C00A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/FC5E9373-B25A-E911-96C8-001E67DDC254.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/FAF4E7DD-4B5A-E911-AECB-1866DAEB40CC.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/FAD3B00E-A159-E911-9022-001E67DDC24A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/FA312F6E-D25A-E911-9170-F01FAFD9C9D0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/E46E4B8A-1E5A-E911-BCE4-AC1F6BAC815A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/F61793D4-7B5A-E911-9C21-0025905B8590.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/F4859CC1-915A-E911-B94C-0242AC130002.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/F2E1B397-D259-E911-B173-842B2B6F5F37.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/F22726F3-9F5A-E911-8FFA-AC1F6BAC8070.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/F2C58245-3F5A-E911-9526-EC0D9A822666.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/F08CF400-9E59-E911-8617-008CFAE452D0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/F07EE92E-795A-E911-A55E-1418774A2299.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/EC4A5628-6F59-E911-8D9C-20CF307C98FC.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/E2C0B956-A65B-E911-908A-A4BF0101DD64.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/E0DAE1AE-915A-E911-9B49-5065F3818271.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/E0AE33DA-A059-E911-8CB5-0242AC130002.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/E09F3227-A45A-E911-A581-0CC47A7C34E6.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/E060E19B-6A5A-E911-91DE-FA163E2C856A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/DE7070EC-075C-E911-BA4F-0CC47A4D765E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/DC40921B-FF59-E911-BD0E-246E96D10990.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/DA821ECC-AE5A-E911-BCF9-1866DA85DFA8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/D8CD661B-9C5A-E911-BFD0-0242AC130002.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/D8966FC0-5359-E911-8C0A-1866DA7F9349.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/D62E044B-985A-E911-A450-0CC47A78A414.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/D24FDAA1-6B5B-E911-8D9C-0242D669068D.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/D0A35D62-C65B-E911-866A-24BE05C38CA1.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/D00617DE-9D5A-E911-9D3E-0CC47A4D76C0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/CCDA9695-3C5A-E911-A3CE-008CFAC91EBC.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/CC4763EF-915D-E911-91A0-AC1F6B1AF224.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/CA42F5C6-B059-E911-99D6-20CF3027A6B0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/C891D540-6459-E911-89A3-0242AC130002.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/C80B56E2-765A-E911-B15B-0CC47A4C8E2E.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/C80598FD-E25A-E911-BCB0-0242AC130002.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/C6DB0B57-285B-E911-A524-FA163E9A9406.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/C6023EE2-2A5B-E911-A124-7CD30AD08CF2.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/C4CD10D8-915D-E911-ADB3-001E67DDBFF7.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/C226ECAD-5659-E911-BFED-0242AC130003.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/C070CCB7-695A-E911-B11A-5065F381C1D1.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/BEEA6F03-495A-E911-8E8E-0242AC130002.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/BE2FADEC-E15A-E911-8081-801844E566F0.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/BCFA9941-375A-E911-AE50-002590D8C7E2.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/BCA5FE61-5C5A-E911-B308-AC1F6BAC7C4A.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/BCA17EA0-D959-E911-A1DA-0242AC130002.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/BC768F4C-475A-E911-86F1-FA163E9C7F3D.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/BC1A7352-D65B-E911-9AA3-549F3525A64C.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/70000/BA744510-925D-E911-A115-AC1F6BABF8D5.root',
],
                  }
    variables = ['recoGenParticles_prunedGenParticles__PAT.obj.m_state.p4Polar_.fCoordinates.fPt', 'recoGenParticles_prunedGenParticles__PAT.obj.m_state.pdgId_', 'recoGenParticles_prunedGenParticles__PAT.obj.m_state.p4Polar_.fCoordinates.fM', 'recoGenParticles_prunedGenParticles__PAT.obj.m_state.p4Polar_.fCoordinates.fEta']

    start_b = time.time()    
    background_list = process_the_samples(backgrounds, None, ['Events'])
    hist_background = get_histograms(background_list, variables, None)

    write_hists_to_file(hist_background, './output_background_b_hists_eta_'+date+'.root') 
    stop_b = time.time()

    signal_list = process_the_samples(signals, 5, ['Events'])
    hist_signal = get_histograms(signal_list, variables, None)

    write_hists_to_file(hist_signal, './output_signal_b_hists_eta_'+date+'.root')  
    stop_s = time.time()

    print "background: ", stop_b - start_b
    print "signal:     ", stop_s - stop_b
    print "total:      ", stop_s - start_b
 
    print 'finished writing'
