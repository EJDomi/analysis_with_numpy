from time import strftime, localtime
import argparse as arg
import os, re
from file_table_functions import *
import ROOT as rt
import numpy as np
import root_numpy as rnp
from collections import OrderedDict
import pickle

date = strftime('%d%b%y', localtime())


def write_sub(output_dir, src_file, file_list, ofile_name, queue_num):
    fsrc = open(os.path.join(output_dir, "submit", src_file+'.submit'), 'w')
    fsrc.write("universe = vanilla \n")
    fsrc.write("executable = " + os.path.join(output_dir, "exe", "run_python_on_condor.sh")+" \n")
    fsrc.write("input = file_table_functions.py\n")
    fsrc.write("getenv = True \n")
    fsrc.write("Arguments = ");
    fsrc.write(os.path.join(output_dir, "file_lists", file_list+"_$(Process).pkl") + " ")
    fsrc.write(os.path.join(output_dir, "out_files", ofile_name)+" \n")
    fsrc.write("output = " + os.path.join(output_dir, "out", ofile_name.replace(".root", ".out")+" \n"))
    fsrc.write("error = " + os.path.join(output_dir, "err", ofile_name.replace(".root", ".err")+" \n"))
    fsrc.write("log = " + os.path.join(output_dir, "log", ofile_name.replace(".root", ".log")+" \n"))
    fsrc.write("Requirements = (Machine != \"red-node000.unl.edu\")\n")
    fsrc.write("queue "+ str(queue_num) + " \n")
    #fsrc.write("cd "+RUN_DIR+" \n")
    #fsrc.write("source ../RestFrames/setup_RestFrames.sh \n")
    fsrc.close()

def write_sh(work_dir, out_dir):
    script_template = '''#!/bin/bash

workdir=WORKDIR
rundir=$(pwd)

cd ${workdir}
eval `scramv1 runtime -sh`
cd ${rundir}

python ${workdir}/make_susy_histograms_p_batch.py --file_list $1 --out_file $2
'''
    out_script = re.sub('WORKDIR', work_dir, script_template)

    with open(os.path.join(out_dir, 'exe', 'run_python_on_condor.sh'), 'w') as f:
        f.write(out_script)



if __name__ == "__main__":
    signals = { 
    'SMS-T2bW_dM' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    #'/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/SMS-T2bW_X05_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_102X.root',
], 
    'SMS-T2-4bd_420' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/erichjs/work/Ewkinos/reducer/CMSSW_10_1_4_patch1/src/KUEWKinoAnalysis_dev_v2/output_10Sep19/Fall17_102X_SMS_I/SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
],
    'SMS-T2-4bd_490' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/erichjs/work/Ewkinos/reducer/CMSSW_10_1_4_patch1/src/KUEWKinoAnalysis_dev_v2/output_10Sep19/Fall17_102X_SMS_I/SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
],
#    'SMS-TChiWH' : [ 
#                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/SMS-TChiWH_WToLNu_HToBB_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_102X',
#],
#    'SMS-T2tt_dM' : [
#                    #'/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/Fall17_102X_SMS_Stop/SMS-T2tt_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/SMS-T2tt_dM-10to80_genHT-160_genMET-80_mWMin-0p1_TuneCP2_13TeV-madgraphMLM-pythia8_Fall17_102X.root',
#], 
              }
    backgrounds = {
    'TTJets_2017' : [
#                     '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X'
#                     '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X'
#                     '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/Fall17_102X_bkgextra/TTTT_TuneCP5_PSweights_13TeV-amcatnlo-pythia8_Fall17_102X',
                     '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                     '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                     '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
],
    'WJets_2017' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
],
    'ZJets_2017' : [
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/ZJetsToNuNu_HT-100To200_13TeV-madgraph_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/ZJetsToNuNu_HT-200To400_13TeV-madgraph_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/ZJetsToNuNu_HT-400To600_13TeV-madgraph_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/ZJetsToNuNu_HT-600To800_13TeV-madgraph_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/ZJetsToNuNu_HT-800To1200_13TeV-madgraph_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph_Fall17_102X',
                    '/home/t3-ku/crogan/NTUPLES/NANO/OLD_29_11_19/NoHadd/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph_Fall17_102X',
],
#    'DY_M50_2017' : [
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/DYJetsToLL_M-50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X',
#],
#    'WW_2017' : [
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WWTo4Q_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WWG_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_102X',
#],
#    'ZZ_2017' : [
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/ZZTo2L2Nu_13TeV_powheg_pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/ZZTo2Q2Nu_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/ZZTo4L_13TeV_powheg_pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_102X',
#],
#    'WZ_2017' : [
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WZG_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WZZ_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_102X',
#                    '/home/t3-ku/crogan/NTUPLES/NANO/NoHadd/Fall17_102X_bkg/WZ_TuneCP5_13TeV-pythia8_Fall17_102X',
#],
                  }

    files_per_job = 100

    background_list = process_the_samples(backgrounds, None, None)
    signal_list = process_the_samples(signals, None, None)

    file_name_base = 'output_{}_risr_95_mixed_'+date+'_{}.root'
  
    working_dir = os.getcwd()

    output_dir = os.path.join(working_dir, 'output_condor_hists_'+date)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
        os.mkdir(os.path.join(output_dir, 'out'))
        os.mkdir(os.path.join(output_dir, 'err'))
        os.mkdir(os.path.join(output_dir, 'log'))
        os.mkdir(os.path.join(output_dir, 'out_files'))
        os.mkdir(os.path.join(output_dir, 'file_lists'))
        os.mkdir(os.path.join(output_dir, 'submit'))
        os.mkdir(os.path.join(output_dir, 'exe'))

    write_sh(working_dir, output_dir)

    list_of_bkg_files = OrderedDict()
    for sample in background_list:
        list_of_bkg_files[sample] = []
        tmp_list = background_list[sample]['files']

        for i in xrange(0, len(tmp_list), files_per_job):
            list_of_bkg_files[sample].append(tmp_list[i:i+files_per_job])

    list_of_sig_files = OrderedDict()
    for sample in signal_list:
        list_of_sig_files[sample] = []
        tmp_list = signal_list[sample]['files']

        for i in xrange(0, len(tmp_list), files_per_job):
            list_of_sig_files[sample].append(tmp_list[i:i+files_per_job])
    
    for sample in list_of_bkg_files:
        submit_list = OrderedDict()
        submit_list[sample] = OrderedDict()
        submit_list[sample]['trees'] = background_list[sample]['trees']
        submit_list[sample]['files'] = None

        src_file = 'hists_{}'.format(sample)
        file_list = '{}_files'.format(sample)

        submit_name = file_name_base.format('background_'+sample, '$(Process)')

        queue_num = len(list_of_bkg_files[sample])

        write_sub(output_dir, src_file, file_list, submit_name, queue_num)

        for i_set, file_set in enumerate(list_of_bkg_files[sample]):
            submit_list[sample]['files'] = file_set

            with open(os.path.join(output_dir, "file_lists", file_list+"_"+str(i_set)+".pkl"), "wb") as f:
                pickle.dump(submit_list, f)
                f.close()

    for sample in list_of_sig_files:
        submit_list = OrderedDict()
        submit_list[sample] = OrderedDict()
        submit_list[sample]['trees'] = signal_list[sample]['trees']
        submit_list[sample]['files'] = None

        src_file = 'hists_{}'.format(sample)
        file_list = '{}_files'.format(sample)

        submit_name = file_name_base.format('signal_'+sample, '$(Process)')

        queue_num = len(list_of_sig_files[sample])

        write_sub(output_dir, src_file, file_list, submit_name, queue_num)

        for i_set, file_set in enumerate(list_of_sig_files[sample]):
            submit_list[sample]['files'] = file_set

            with open(os.path.join(output_dir, "file_lists", file_list+"_"+str(i_set)+".pkl"), "wb") as f:
                pickle.dump(submit_list, f)
                f.close()
        

 
