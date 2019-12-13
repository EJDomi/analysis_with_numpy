#!/usr/bin/env python

"""
New thing to try out, doing analysis with numpy, converting back to ROOT to do histogramming
Creator: Erich Schmitz
Date: Feb 22, 2019
"""

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
            
        for ifile, in_file in enumerate(list_of_files_[sample]['files']):
            for tree_name in list_of_files_[sample]['trees']:
                sample_array = get_tree_info_singular(sample, in_file, tree_name, variable_list_, cuts_to_apply_)
                if sample_array is None: continue

                print '\nGetting Histograms for:', sample, tree_name, in_file
                print 'file: ', ifile+1, ' / ', len(list_of_files_[sample]['files'])

                gen_pdgid = np.array(sample_array['GenPart_pdgId'])
                gen_pt = np.array(sample_array['GenPart_pt'])
                gen_mass = np.array(sample_array['GenPart_mass'])

                is_susy = np.full(len(gen_pt), True)

                if 'SMS_T2bW' in sample:
                    stop_mass = np.array([mass[np.logical_or(pdgid == 1000006, pdgid == 2000006)] for mass, pdgid in zip(gen_mass, gen_pdgid)])
                    child_mass = np.array([mass[pdgid == 1000022] for mass, pdgid in zip(gen_mass, gen_pdgid)])
                    is_susy = np.array([True if np.any(stop == 500.) and np.any(child == 460.) else False for stop, child in zip(stop_mass, child_mass)])
                
 
                b_pt = np.array([pt[np.abs(pdgid) == 5] for pt, pdgid, susy in zip(gen_pt, gen_pdgid, is_susy) if susy])

                if not np.any(is_susy): continue

                rnp.fill_hist(hist[sample][tree_name]['PT_b'], np.concatenate(b_pt))

                print 'finished filling'
    return hist


             

if __name__ == "__main__":

    signals = { 
    'SMS_T2bW_460' : [
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/F6E8C114-304C-6F4A-B7B6-AF9D593F0DCA.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/F39E109A-4BDE-6145-99C6-BCE759524078.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/ECB66637-A683-7648-8C3D-1A2CDA50E34C.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/E55091C9-41B6-FE43-836A-B8C530B2F1F3.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/CED699B1-00E1-4640-A305-5D716028D980.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/CE93E383-DD0D-1E44-94D9-F47BE3AA10E2.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/CD05AD01-60F0-274F-84E9-8B9894F321CD.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/C53A7A7A-CDA4-E649-BDFD-599C700D2FDF.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/BAA5B2A1-D930-C04E-A05C-46805BE815BD.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/AEBFA8FA-50F9-D34A-8013-6ABF434DC64A.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/9DE5C052-0FCA-8F45-BB01-086382051BBD.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/93B10C45-436A-B449-88C3-2BBD8C36BBBB.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/931E57C7-A5A7-124D-B09A-7A7D98C5CA5E.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/8F778C79-CC04-4D41-8ECE-1D7E34A3C4E4.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/8997097F-1F71-394B-BD6D-2C793B14F080.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/841C9FED-A678-4742-A746-68D4C91FEB8B.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/8211BF93-EF1B-5142-AB5F-9EF6794E1212.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/7E145827-E92D-8145-826D-8B968C7B0668.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/7B129057-525A-0B44-ACFE-17BDC312BD84.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/77C7C522-1DAF-A84A-9007-32DFEE7CECEE.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/4B007784-D494-6A4F-B6A3-0940C227B338.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/40049FEA-AC30-B74E-BF41-47DBED5A432C.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/3E86B736-102F-2A4B-90A6-6DBBA97E25AF.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/3C0CBDFE-DC78-8446-9C22-5F1E0BB74DEA.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/26825472-58FB-0945-958A-DDD5EB4F1F80.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/22D5C15E-464D-FC41-8F64-59EF11D6AE03.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/1ED85FE6-0412-3642-98E7-AE013E34671E.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/0F41DB6F-DA06-1540-AD80-9273B3360DE2.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/0AB1C506-4659-1C44-ABCE-6C431D350920.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/270000/036A9E3E-E927-6445-BE09-4FE0335DFCBD.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/260000/E932AE35-3B2A-7146-9D6D-DC0CA5424C4A.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/260000/E47C1F38-98C7-474B-B342-92ED36DFB16D.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/260000/D0FCACC6-E744-2F46-AA5D-A8E7C55D34DE.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/260000/BEC623AA-9913-014A-8CD9-1BD0B26DE8AD.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/260000/ADC9F442-82B3-DB46-9F12-5BC824A20B59.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/260000/6225B7D0-8D07-914C-8DEC-F6BE3CF277BE.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/260000/3DDBD4FD-E0B8-C443-9067-0EE3A97B3EDE.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/260000/2BD7E015-5503-C34C-B3ED-1528FA0ECF13.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/260000/11DE276C-C724-6F43-B3F7-C8C519BE9241.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/260000/0A6250E0-4683-994B-AB03-CC3D98AB5F83.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/FA56C5D5-2DE7-2846-8E37-11C1F37B0DFA.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/DFBAFA01-D653-9747-AE8C-53DCB8A97AF8.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/DE645841-2574-5840-A0C7-6942A6E6E940.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/B600ED3B-70BB-0A48-BD57-EDE9A3A9AD7B.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/A4F1AF57-A218-444D-87BB-D6944A89BF0F.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/9DD9E8C0-2E44-1F4C-AA33-905A1CF6CBBC.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/95CD1334-44C3-6D44-9962-5D3CDA719C1B.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/8FD26675-2F17-C741-A988-A29054F9FB30.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/8EBEF0A9-A5AB-E44C-B6EE-4DF479765B78.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/8D4E5174-65EA-6340-9FBA-5DAF20E18501.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/87FC3F0E-830E-2A40-9C2B-3F5FCCBFE755.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/7FBF2FE6-20F0-7E41-BE17-DF07AB6C95C7.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/7AD28964-4690-8A48-88AB-473FD543269E.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/72124704-E0A8-9946-9B5D-F577AAD3CEE3.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/63D03595-FE14-A24D-8042-164B32BDA531.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/38482234-CEDE-404E-8AF5-26736C72CD8A.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/25AB23C8-D568-5A49-9EB1-D60146AC8CFE.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/24FFEE75-8382-DD4D-A0B7-BA4443AE2763.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/17931266-34EC-C749-B939-58FECFDA3AAD.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/10DDA39D-D6F0-6A47-BFF6-3252461EFFE4.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/FA60944C-BE5F-FC41-AC0E-6024A6E7C7D5.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/FA06E916-C237-8741-A244-2A3483C620C9.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/EFB5EE75-9A19-3142-9DF8-26FA040CE30A.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/EC7D4070-FF72-4946-A3E8-50BF9AC17929.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/E5F6C47F-B593-E049-9AE8-04F4A6B29D32.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/E0A2BDE9-8C84-3448-8B04-7F39331FF9E1.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/D9B1A1AA-6144-8B4C-A397-8DA40A3D710E.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/CB26EB79-3F76-CF44-9DDA-158C7DC5F3DB.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/CA7D457D-D9F2-384D-879E-5FB9EC7DCBC5.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/BA00ABC5-19AA-CD4D-A67C-1C7E953680B3.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/9D5132F3-C3D7-4045-9536-25EE641133FD.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/8FEECC10-9C5E-1E41-A9E6-3801919C6249.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/8FADA7D7-AF46-1D4F-B25D-0CF4C8D76CE1.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/72BF2D6A-BE37-774B-B658-C447B8CB754D.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/6C6B2725-6CDD-BC48-909A-6FFF924AE180.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/6009900F-D7E4-CF4F-934A-DF5FBA8F0E81.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/584AC408-3AEA-694F-B316-5D419A49080C.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/36AC4831-B642-2441-B3C7-664CF4AD0B20.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/35E7EE77-A3A4-8343-9B54-2167E7E7EE78.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/2F6EC6FC-FAF2-AF48-A5E6-559541F76BEB.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/25A5FE10-E5E7-C44E-B827-264AE2CB9BA1.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/1E08AB65-4FEE-CE45-A2B2-761A9B59A39A.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/1DFDD190-1F39-834E-8B08-175B7B5B48CC.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/171FD31B-7794-9F49-8562-6B9505EDA71F.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUFall17Fast_Nano1June2019_102X_mc2017_realistic_v7-v1/100000/0691CA8F-361A-984C-BDA3-D051794163D5.root',
],
    'SMS_T2_4bd_420' : [
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/CBD07CFB-38D0-124C-A1B2-140586FCA86F.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP2_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/6712BCF1-402B-3140-8F1F-68FA8281A9FD.root',
],
    'SMS_T2_4bd_490' : [
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/EC9411D3-C6DB-0C48-AB2E-5F7DD0F42B4B.root',
'root://cmsxrootd.fnal.gov////store/mc/RunIIFall17NanoAODv5/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/7D503526-1D50-0C47-8C14-2EB9882BAC65.root',
],
              }
    backgrounds = {
    'TTJets' : [
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17NanoAODv5/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/8D625D30-C72D-0747-B45F-FDDFC7ABE4A2.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17NanoAODv5/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/23AE3FB4-D814-E544-85B8-7B001BBA8416.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17NanoAODv5/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/9A5C1FCA-43DD-B244-B58A-9D277A8A4698.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17NanoAODv5/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/687DB9BB-82D1-2C49-B967-B2508D237ADA.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17NanoAODv5/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/160BBBD6-645E-FE4E-BB4C-4A657C8BCD11.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17NanoAODv5/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/6844D099-CC6F-7940-8455-23185E0308C8.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17NanoAODv5/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/F902413C-35DF-624C-9ECC-C5B4F8BDBF08.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17NanoAODv5/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/80455F65-985C-AA4D-A6FC-0EB1B55C2303.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17NanoAODv5/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/C541F35D-DB04-C448-8687-75737F0FF800.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17NanoAODv5/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/B4642F0E-7FFC-1D46-93E3-11B59C616558.root',
'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17NanoAODv5/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/09E97185-D9B7-5846-BC6F-7D4DB0B5403E.root',
],
    'WJets_2017' : [
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/2C8EC1C3-0127-DF4D-AFD3-2EA0C7D89BC5.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/0B2720E3-5659-4E4F-B79A-FE3F778DD587.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/837DEDDE-9F31-8F41-A47F-F1DADDDE9764.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/39C304D6-11BE-4B4D-938F-445B8D5DB162.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/5EFBD5BC-FC04-454F-9347-FF6F152F5E03.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/8C7C3BA1-0878-A14F-93A5-9ACF79A5D0D7.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/E2BBA38E-15CF-594E-B5CF-E1B0A686F11E.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/D6FBBA28-8242-9841-864B-66AB96E27953.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/B9095E2A-5E71-E945-B510-38A005C8CE9C.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/120000/F6E7160B-28F6-9343-ACD2-3DCF8CBD141A.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/42CC51F7-8EA2-034C-A576-7528FEFC440D.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/6F1C2DAA-252C-D04F-8840-EDF02C313DA1.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/4DA97933-536F-284E-873B-FCEC4BA49EF0.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/5FF7592A-F419-4047-9D8C-92ECA7482194.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/DB82F5DF-6B2D-2A4D-834E-B790B72849B3.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/DDE76EE1-BCE9-5249-82AD-90886AB7811E.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/898A0E77-53B8-C047-BD94-F0C66768A9F8.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/BCF14775-8055-0849-B6F3-64F6F4B3259E.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/59A3A9A4-660B-7B4F-B858-23CFE7C09303.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/85CAE966-4A71-A249-9AFA-92E8970B9795.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/34D47303-0A92-074E-9C5B-AD971708ED96.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/110000/941CE317-60F7-0E42-AAD4-D56DAE581325.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/100000/757BEDDA-9CCB-6948-9B6A-A6E3CD888DF6.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/70000/A59212FC-B70E-CC43-B7AD-8854599F05FC.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/70000/CABE3F37-E573-6B4B-9B25-08A65E397CD0.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/70000/12FFEDDF-9119-2A44-B797-1BFA8CB3541E.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/70000/D32C8308-BA33-EA4C-B567-F417369C1F01.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/70000/865FEC14-C764-B245-B1DB-02246F3C27C8.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/70000/146B12A0-FDC4-C545-816F-21619D0B554F.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/70000/6540338A-0F5B-6540-BF68-F13AE1B782BE.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/A688C860-9E5C-D540-88F7-F76678F8B427.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/46953206-988B-A449-B90B-272646857B93.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/100000/058842A0-738B-974B-BD82-2AF46391D30D.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/70000/098CF623-BACD-1E4F-BCED-71ABBFD5AB71.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/70000/3C291788-7B0F-1D48-AD7D-1615370C2C71.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/70000/9FE8AB8D-997D-E04D-9B1B-72C1698EB053.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/70000/B58F2875-E116-8A43-8AD0-AFCDB5D4712B.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/110000/314170BC-2FEB-BF4C-BD84-600066E65515.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/110000/38826084-46D1-974E-8B30-1035BDBCE3D8.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/110000/674D5AA5-AF4A-9B4C-94A4-188B09623820.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/AEF27C79-CC77-CF4F-9E0B-8331A2F785CC.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/AA820153-7239-C740-9522-69D750E8D607.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/041F214F-53B9-0D45-8941-2F2DD4A3F437.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/4ED071AB-5AF5-684F-AAC6-33816E69CB80.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/D04BEE48-672E-9D4A-8693-B1C25E6590AF.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/32747F8B-9047-E242-A9ED-3B29B934D985.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/D44F6647-C06C-2649-B573-7E51112F280F.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/B668E5F4-CBC7-B048-AEAB-F5AEFAC05797.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/18F7AB4F-E8F0-1E43-956B-AD2DEE1A3D60.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/D0C804C2-BF94-7744-A61C-046659DD63AA.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/E205F637-0B8B-AD40-80FC-379F90703ECE.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/B84A16DE-4C1A-8148-9780-6DCAD270A90A.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/99C73928-3CAA-354C-A6D8-5618188D3D3D.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/60000/7CFA0519-E5C7-4943-BEB4-929EED2CCD11.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/130000/44FF5CB9-63CB-104A-B7D1-8F1D50B43967.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/130000/CCDC79A5-A457-C745-8726-AB673462A52A.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/130000/6EF3EE57-131E-454E-B774-BBAE54F816A2.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/F5DA139A-A45B-1648-A3AC-DAFFFA0A5B4F.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/463398E0-D65A-1F45-A379-47BBF3770646.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/EF5FA0BC-3978-3D49-B164-D49AF8D8CC72.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/64B44C95-57FE-594C-AE56-A56B909B0A31.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/3623656D-16F7-CA43-B49F-B52B33A193B8.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/07935E1A-5F24-8243-96C8-15B66D9957FB.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/5D6BD219-1FDD-C247-BDF2-F0603999B325.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/43FE9772-5EC9-C34D-B362-89009BE3BD56.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/B3DD15EE-FA9B-0545-AD6B-A185A2DA68DC.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/11BA0BCA-327E-4648-BDA3-A002EC799256.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/E2D84EDD-C95E-BE44-B417-791CBF53E3A6.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/2E736E6C-56C7-6641-B5BC-02421492F165.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/E783C8A4-4C29-4043-9A38-47AE52879B52.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/AA46E5AF-6251-0442-BB0A-F49B8ED69A2D.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/DD2F1FC6-90DD-8341-984E-2D8F1B79A04F.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/260000/11A6A3EB-2DE9-7444-A98D-F9CD2AD83C15.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/130000/2D93BA25-A2E8-ED40-817A-11A24BF986B8.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_new_pmx_102X_mc2017_realistic_v7-v1/130000/F528B56A-228A-C647-9013-B4943E0F7B3B.root',
],
                  }
    variables = ['GenPart_pdgId', 'GenPart_pt', 'GenPart_mass']

    start_b = time.time()    
    #background_list = process_the_samples(backgrounds, None, ['Events'])
    #hist_background = get_histograms(background_list, variables, None)

    #write_hists_to_file(hist_background, './output_background_b_hists.root') 
    stop_b = time.time()

    signal_list = process_the_samples(signals, None, ['Events'])
    hist_signal = get_histograms(signal_list, variables, None)

    write_hists_to_file(hist_signal, './output_signal_b_hists.root')  
    stop_s = time.time()

    print "background: ", stop_b - start_b
    print "signal:     ", stop_s - stop_b
    print "total:      ", stop_s - start_b
 
    print 'finished writing'
