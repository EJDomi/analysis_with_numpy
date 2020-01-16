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
import pickle
import argparse as arg

rt.gROOT.SetBatch()
rt.TH1.AddDirectory(rt.kFALSE)

date = strftime('%d%b%y', localtime())
########## get_histograms function template ############

########################################################


def make_tables(list_of_files_, variable_list_, cuts_to_apply_=None):
    
    hist = OrderedDict()
    reference = [
            # s-style regions
            'total',
            '2l',
            '2l_ge1st',
            '2l_ge1sv',
            '1l',
            '1l_0ssv',
            '1l_0sb',
            '1l_1sb',
            '1l_ge1sb',
            '1l_ge1v',
            '1l_ge1st',
            '1l_ge1sj',
            'ge1sv',
            'ge2sv',
            'ge1sb',

            '1_el',
            '1_el_ge1st',
            '1_el_0isrb',
            '1_el_0isrb_ge1sj',
            '1_el_0isrb_ge1st',
            '1_mu',
            '1_mu_ge1st',
            '1_mu_0isrb',
            '1_mu_0isrb_ge1sj',
            '1_mu_0isrb_ge1st',
            # isr-style regions

            '1l_0isrb_ge1sb_ge1sv',
            '1l_0isrb_ge1sv_philmetg',
            '1l_0isrb_ge1sv_philmetl',

            '1l_0isrb_ge1sv_mscut1_philmetl',
            '1l_0isrb_ge1sv_mscut1_philmetg',

            '1l_0isrb_ge1sv_mscut3_philmetg',

            '1l_0isrb_ge1st_mscut1',
            '1l_0isrb_ge1st_mscut2',
            '1l_0isrb_ge1st_mscut3',

            '1l_0isrb_ge1sb_mscut1',
            '1l_0isrb_ge1sb_mscut2',
            '1l_0isrb_ge1sb_mscut3',

            '1l_0isrb_ge1sv_mscut1',
            '1l_0isrb_ge1sv_mscut2',
            '1l_0isrb_ge1sv_mscut3',

            '1l_0isrb_ge1v_mscut1',

            '1l_0isrb_ge2sv',
            '1l_0isrb_ge2st',

            '1l_0isrb_ge1st',
            '1l_0isrb_ge1st_philmetg',
            '1l_0isrb_ge1st_philmetl',

            'isr_ge1v',
            'isr_ge2v',
            'isr_ge1b',
            '1l_1isrb',
            '1l_ge1isrb',
            '1l_ge1isrt',
            '1l_ge2isrb',
            '2l_0isrb',
            '2l_ge1isrb',
            '2l_0isrb_ge1st',
            '2l_0isrb_ge1sv',
            '2l_0isrb_r8',
            '2l_0isrb_r9',
            '2l_0isrb_r95',
            '2l_0isrb_ge1st_r8',
            '2l_0isrb_ge1st_r9',
            '2l_0isrb_ge1st_r95',
            '2l_0isrb_ge1sv_r8',
            '2l_0isrb_ge1sv_r9',
            '2l_0isrb_ge1sv_r95',
            '2l_ge1isrb_r8',
            '2l_ge1isrb_r9',
            '2l_ge1isrb_r95',
            '1l_0isrb',
            '1l_0isrt',
            '1l_0isrb_0ssv',
            '1l_0isrb_0s_r8',
            '1l_0isrb_0s_r9',
            '1l_0isrb_0s_r95',
            '1l_ge2isrb_r8',
            '1l_ge2isrb_r9',
            '1l_ge2isrb_r95',
            '1l_0isrb_ge1sj_m1_r95',
            '1l_0isrb_ge1sj',
            '1l_0isrb_ge1sb',
            '1l_0isrb_ge1sv',
            '1l_0isrb_ge1v',
            '1l_0isrb_ge1sb_ge1sv_r8',
            '1l_0isrb_ge1sb_ge1sv_r9',
            '1l_0isrb_ge1sb_ge1sv_r95',
            '1l_0isrb_ge1sj_mscut1',
            '1l_0isrb_ge1sj_mscut2',
            '1l_0isrb_ge1st_m1_r9',
            '1l_0isrb_ge1st_m1_r95',
            '1l_0isrb_ge1st_m3_r8',
            '1l_0isrb_ge1sv_m1_r95',
            '1l_0isrb_ge1sv_m1_r9',
            '1l_0isrb_ge1sv_mscut3',
            '1l_0isrb_ge1sv_m1_philmetl_r95',
            '1l_0isrb_ge1sv_m3_philmetg_r8',
            '1l_0isrb_ge2st_r8',
            '1l_0isrb_ge2st_r9',
            '1l_0isrb_ge2st_r95',
            '1l_0isrb_ge1sb_m1_r9',
            '1l_0isrb_ge1sb_m1_r95',
            '1l_0isrb_ge1sb_m3_r8',
            '1l_0isrb_ge1sv_pt12_r8',
            '1l_0isrb_ge1sv_pt12_r9',
            '1l_0isrb_ge1sv_pt12_r95',
            '1l_0isrb_ge1sv_pt10_r8',
            '1l_0isrb_ge1sv_pt10_r9',
            '1l_0isrb_ge1sv_pt10_r95',
            '1l_0isrb_ge1sv_pt7_r8',
            '1l_0isrb_ge1sv_pt7_r9',
            '1l_0isrb_ge1sv_pt7_r95',
            '1l_0isrb_ge1sv_dphlep_r8',
            '1l_0isrb_ge1sv_dphlep_r9',
            '1l_0isrb_ge1sv_dphlep_r95',
            '1l_0isrb_ge1sv_dphlep_pt7_r8',
            '1l_0isrb_ge1sv_dphlep_pt7_r9',
            '1l_0isrb_ge1sv_dphlep_pt7_r95',
            '1l_0isrb_ge1sv_chp_r8',
            '1l_0isrb_ge1sv_chp_r9',
            '1l_0isrb_ge1sv_chp_r95',

            # weighted entries
            # s-style regions
            'total_w',
            '2l_w',
            '2l_ge1st_w',
            '2l_ge1sv_w',
            '2l_0isrb_w',
            '2l_ge1isrb_w',
            '2l_0isrb_ge1st_w',
            '2l_0isrb_ge1sv_w',
            '2l_0isrb_r8_w',
            '2l_0isrb_r9_w',
            '2l_0isrb_r95_w',
            '2l_0isrb_ge1st_r8_w',
            '2l_0isrb_ge1st_r9_w',
            '2l_0isrb_ge1st_r95_w',
            '2l_0isrb_ge1sv_r8_w',
            '2l_0isrb_ge1sv_r9_w',
            '2l_0isrb_ge1sv_r95_w',
            '2l_ge1isrb_r8_w',
            '2l_ge1isrb_r9_w',
            '2l_ge1isrb_r95_w',

            '1l_w',
            '1l_0ssv_w',
            '1l_0sb_w',
            '1l_1sb_w',
            '1l_ge1sb_w',
            '1l_ge1v_w',
            '1l_ge1st_w',
            '1l_ge1sj_w',
            'ge1sv_w',
            'ge2sv_w',
            'ge1sb_w',

            'isr_ge1v_w',
            'isr_ge2v_w',
            'isr_ge1b_w',
            '1_el_w',
            '1_el_ge1st_w',
            '1_el_0isrb_w',
            '1_el_0isrb_ge1sj_w',
            '1_el_0isrb_ge1st_w',
            '1_mu_w',
            '1_mu_ge1st_w',
            '1_mu_0isrb_w',
            '1_mu_0isrb_ge1sj_w',
            '1_mu_0isrb_ge1st_w',

            '1l_0isrb_ge1sb_ge1sv_w',
            '1l_0isrb_ge1sv_philmetg_w',
            '1l_0isrb_ge1sv_philmetl_w',

            '1l_0isrb_ge1sv_mscut1_philmetl_w',
            '1l_0isrb_ge1sv_mscut1_philmetg_w',

            '1l_0isrb_ge1sv_mscut3_philmetg_w',

            '1l_0isrb_ge1st_mscut1_w',
            '1l_0isrb_ge1st_mscut2_w',
            '1l_0isrb_ge1st_mscut3_w',

            '1l_0isrb_ge1sb_mscut1_w',
            '1l_0isrb_ge1sb_mscut2_w',
            '1l_0isrb_ge1sb_mscut3_w',

            '1l_0isrb_ge1sv_mscut1_w',
            '1l_0isrb_ge1sv_mscut2_w',
            '1l_0isrb_ge1sv_mscut3_w',

            '1l_0isrb_ge1v_mscut1_w',

            '1l_0isrb_ge2sv_w',
            '1l_0isrb_ge2st_w',

            '1l_0isrb_ge1st_w',
            '1l_0isrb_ge1st_philmetg_w',
            '1l_0isrb_ge1st_philmetl_w',
            '1l_1isrb_w',
            '1l_ge1isrb_w',
            '1l_ge1isrt_w',
            '1l_ge2isrb_w',
            '1l_0isrb_w',
            '1l_0isrt_w',
            '1l_0isrb_0ssv_w',
            '1l_0isrb_0s_r8_w',
            '1l_0isrb_0s_r9_w',
            '1l_0isrb_0s_r95_w',
            '1l_ge2isrb_r8_w',
            '1l_ge2isrb_r9_w',
            '1l_ge2isrb_r95_w',
            '1l_0isrb_ge1sj_m1_r95_w',
            '1l_0isrb_ge1sj_w',
            '1l_0isrb_ge1sb_w',
            '1l_0isrb_ge1sv_w',
            '1l_0isrb_ge1v_w',
            '1l_0isrb_ge1sb_ge1sv_r8_w',
            '1l_0isrb_ge1sb_ge1sv_r9_w',
            '1l_0isrb_ge1sb_ge1sv_r95_w',
            '1l_0isrb_ge1sj_mscut1_w',
            '1l_0isrb_ge1sj_mscut2_w',
            '1l_0isrb_ge1st_m1_r9_w',
            '1l_0isrb_ge1st_m1_r95_w',
            '1l_0isrb_ge1st_m3_r8_w',
            '1l_0isrb_ge1sv_m1_r95_w',
            '1l_0isrb_ge1sv_m1_r9_w',
            '1l_0isrb_ge1sv_mscut3_w',
            '1l_0isrb_ge1sv_m1_philmetl_r95_w',
            '1l_0isrb_ge1sv_m3_philmetg_r8_w',
            '1l_0isrb_ge2st_r8_w',
            '1l_0isrb_ge2st_r9_w',
            '1l_0isrb_ge2st_r95_w',
            '1l_0isrb_ge1sb_m1_r9_w',
            '1l_0isrb_ge1sb_m1_r95_w',
            '1l_0isrb_ge1sb_m3_r8_w',
            '1l_0isrb_ge1sv_pt12_r8_w',
            '1l_0isrb_ge1sv_pt12_r9_w',
            '1l_0isrb_ge1sv_pt12_r95_w',
            '1l_0isrb_ge1sv_pt10_r8_w',
            '1l_0isrb_ge1sv_pt10_r9_w',
            '1l_0isrb_ge1sv_pt10_r95_w',
            '1l_0isrb_ge1sv_pt7_r8_w',
            '1l_0isrb_ge1sv_pt7_r9_w',
            '1l_0isrb_ge1sv_pt7_r95_w',
            '1l_0isrb_ge1sv_dphlep_r8_w',
            '1l_0isrb_ge1sv_dphlep_r9_w',
            '1l_0isrb_ge1sv_dphlep_r95_w',
            '1l_0isrb_ge1sv_dphlep_pt7_r8_w',
            '1l_0isrb_ge1sv_dphlep_pt7_r9_w',
            '1l_0isrb_ge1sv_dphlep_pt7_r95_w',
            '1l_0isrb_ge1sv_chp_r8_w',
            '1l_0isrb_ge1sv_chp_r9_w',
            '1l_0isrb_ge1sv_chp_r95_w'
    ]

    for sample in list_of_files_:
        hist[sample] = OrderedDict()
        for tree_name in list_of_files_[sample]['trees']:
            print '\nReserving Tables for:', sample, tree_name 
            # Reserve tables

            n_regions = len(reference)
            hist[sample][tree_name] = np.zeros(n_regions, dtype=np.float)

        for ifile, in_file in enumerate(list_of_files_[sample]['files']):
            for tree_name in list_of_files_[sample]['trees']:
                sample_array = get_tree_info_singular(sample, in_file, tree_name, variable_list_, cuts_to_apply_)
                if sample_array is None: continue

                print '\nGetting Histograms for:', sample, tree_name, in_file
                print 'file: ', ifile+1, ' / ', len(list_of_files_[sample]['files'])

                risr = np.array(sample_array['RISR'])
                if len(risr) ==0: continue
                ptisr = np.array(sample_array['PTISR'])
                ms = np.array(sample_array['MS'])
                nlep = np.array(sample_array['Nlep'])
                nlep_s = np.array(sample_array['Nlep_S'])
                nlep_isr = np.array(sample_array['Nlep_ISR'])
                njet_s = np.array(sample_array['Njet_S'])
                njet_isr = np.array(sample_array['Njet_ISR'])
                nbjet_s = np.array(sample_array['Nbjet_S'])
                nbjet_isr = np.array(sample_array['Nbjet_ISR'])
                nsv_s = np.array(sample_array['NSV_S'])
                nsv_isr = np.array(sample_array['NSV_ISR'])

#                pt_jet = np.array(sample_array['PT_jet'])
#                pt_sv = np.array(sample_array['PT_SV'])
#                isr_index_jet = np.array(sample_array['index_jet_ISR'])
#                s_index_jet = np.array(sample_array['index_jet_S'])
                isr_index_sv = np.array(sample_array['index_SV_ISR'])
                s_index_sv = np.array(sample_array['index_SV_S'])
#                bjet_tag = np.array(sample_array['Btag_jet'])

                pt_lep = np.array(sample_array['PT_lep'])
                pa_lep = np.array(sample_array['PLa'])
                pb_lep = np.array(sample_array['PLb'])
                ch_lep = np.array(sample_array['Charge_lep'])
                phi_lep = np.array(sample_array['Phi_lep'])
                mini_lep = np.array(sample_array['MiniIso_lep'])
                sip3d_lep = np.array(sample_array['SIP3D_lep'])
                id_lep = np.array(sample_array['ID_lep'])
                pdgid_lep = np.array(sample_array['PDGID_lep'])
                dphi_lep_s = np.array(sample_array['dphi_lep_S'])

                eta_sv = np.array(sample_array['Eta_SV'])

                risr = np.array([entry[:2] for entry in risr])
                ptisr = np.array([entry[:2] for entry in ptisr])
                ms = np.array([entry[:2] for entry in ms])
                njet_s = np.array([entry[:2] for entry in njet_s])
                njet_isr = np.array([entry[:2] for entry in njet_isr])
                nbjet_s = np.array([entry[:2] for entry in nbjet_s])
                nbjet_isr = np.array([entry[:2] for entry in nbjet_isr])
                nsv_s = np.array([entry[:2] for entry in nsv_s])
                nsv_isr = np.array([entry[:2] for entry in nsv_isr])
                #ptcm = np.array([entry[:2] for entry in ptcm])
                #dphi = np.array([entry[:2] for entry in dphi])
#                isr_index_jet = np.array([entry[:2] for entry in isr_index_jet])
#                s_index_jet = np.array([entry[:2] for entry in s_index_jet])
                isr_index_sv = np.array([entry[:2] if np.shape(entry[:2]) != (0,) else [np.array([],dtype=np.int32),np.array([],dtype=np.int32)] for entry in isr_index_sv])
                s_index_sv = np.array([entry[:2] if np.shape(entry[:2]) != (0,) else [np.array([],dtype=np.int32),np.array([],dtype=np.int32)] for entry in s_index_sv])
                pa_lep = np.array([entry[:2] for entry in pa_lep])
                pb_lep = np.array([entry[:2] for entry in pb_lep])
                dphi_lep_s = np.array([entry[:2] for entry in dphi_lep_s])

                risr = risr[:, 1]
                ptisr = ptisr[:, 1]
                ms = ms[:, 1]
                njet_s = njet_s[:, 1]
                njet_isr = njet_isr[:, 1]
                nbjet_s = nbjet_s[:, 1]
                nbjet_isr = nbjet_isr[:, 1]
                nsv_s = nsv_s[:, 1]
                nsv_isr = nsv_isr[:, 1]
                pa_lep = pa_lep[:, 1]
                pb_lep = pb_lep[:, 1]
                dphi_lep_s = dphi_lep_s[:, 1]

                #ptcm = ptcm[:, 2]
                #dphi = dphi[:, 2]
#                isr_index_jet = isr_index_jet[:, 1]
#                s_index_jet = s_index_jet[:, 1]
                isr_index_sv = isr_index_sv[:, 1]
                s_index_sv = s_index_sv[:, 1]
                #isr_index_lep = isr_index_lep[:, 2]
                #s_index_lep = s_index_lep[:, 2]

                #risr = risr[:, 1]
                #ptisr = ptisr[:, 1]
                #ptcm = ptcm[:, 1]
                #isr_index_jet = isr_index_jet[:, 1]
                #s_index_jet = s_index_jet[:, 1]
                #isr_index_lep = isr_index_lep[:, 1]
                #s_index_lep = s_index_lep[:, 1]
                #dphi = dphi[:, 1]

                # risr_lepV_jetI = risr[:,0]
                # risr_lepV_jetA = risr[:,1]
                # risr_lepA_jetA = risr[:,2]

                met = np.array(sample_array['MET'])
                met_phi = np.array(sample_array['MET_phi'])
                weight = np.array(sample_array['weight'])
                weight = 137. * weight

                eta_sv_s = np.array([eta[ind] for eta, ind in zip(eta_sv, s_index_sv)])
                eta_sv_isr = np.array([eta[ind] for eta, ind in zip(eta_sv, isr_index_sv)])
               
                nsv_eta = np.array([len(eta[np.abs(eta)<1.5]) for eta in eta_sv])
                nsv_eta_s = np.array([len(eta[np.abs(eta)<1.5]) for eta in eta_sv_s])
                nsv_eta_isr = np.array([len(eta[np.abs(eta)<1.5]) for eta in eta_sv_isr])

                p_lep = pa_lep + pb_lep
                one_lep = nlep == 1
                two_lep = nlep == 2
                ################  Choosing Lepton ID #####################
                ################        Medium       #####################
                mini_lep = np.array([mini[lid>=3]*pt[lid>=3] for mini, lid, pt in zip(mini_lep, id_lep, pt_lep)])
                sip3d_lep = np.array([ips[lid>=3] for ips, lid in zip(sip3d_lep, id_lep)])
                pt_lep = np.array([pt[lid>=3] for pt, lid in zip(pt_lep, id_lep)])
                ch_lep = np.array([ch[lid>=3] for ch, lid in zip(ch_lep, id_lep)])
                phi_lep = np.array([phi[lid>=3] for phi, lid in zip(phi_lep, id_lep)])
                pdgid_lep = np.array([pt[lid>=3] for pt, lid in zip(pdgid_lep, id_lep)])
                ################        Tight        #####################
                #pt_lep = np.array([pt[lid>=4] for pt, lid in zip(pt_lep, id_lep)])
                #mini_lep = np.array([mini[lid>=4] for mini, lid in zip(mini_lep, id_lep)])
                ##########################################################
 
                is_mu = np.array([ np.abs(pdg) == 13 for pdg in pdgid_lep])
                is_el = np.array([ np.abs(pdg) == 11 for pdg in pdgid_lep])
                
                pt_mu = np.array([ pt[mu] for pt, mu in zip(pt_lep, is_mu)])
                pt_el = np.array([ pt[el] for pt, el in zip(pt_lep, is_el)])

                mini_mu = np.array([ mini[mu] for mini, mu in zip(mini_lep, is_mu)])
                mini_el = np.array([ mini[el] for mini, el in zip(mini_lep, is_el)])

                sip3d_mu = np.array([ sip3d[mu] for sip3d, mu in zip(sip3d_lep, is_mu)])
                sip3d_el = np.array([ sip3d[el] for sip3d, el in zip(sip3d_lep, is_el)])

#                pt_35_lep = pt_lep
                pt_35_lep = np.array([pt[np.logical_and(ips<4., mini<6.0)] for pt, ips, mini in zip(pt_lep, sip3d_lep, mini_lep)])
#                pt_35_lep = np.array([ pt[np.logical_and(pt>5, mini<0.1)] for pt, mini in zip(pt_lep, mini_lep)])
#                pt_35_mu = pt_mu
#                pt_35_mu = np.array([ pt[mini<6.] for pt, mini in zip(pt_mu, mini_mu)])
                pt_35_mu = np.array([ pt[np.logical_and(ips<4., mini<6.0)] for pt, ips, mini in zip(pt_mu, sip3d_mu, mini_mu)])
#                pt_5_el = pt_el
                pt_5_el = np.array([ pt[np.logical_and(ips<4., mini<6.0)] for pt, ips, mini in zip(pt_el, sip3d_el, mini_el)])
              
#                pt_less_12_lep = np.array([ pt[np.logical_and(pt<12, mini<6.0)] for pt, mini in zip(pt_lep, mini_lep)])

                ch_35_lep = np.array([pt[np.logical_and(ips<4., mini<6.0)] for pt, ips, mini in zip(ch_lep, sip3d_lep, mini_lep)])


                phi_mini_lep = np.array([ phi[np.logical_and(ips<4., mini<6.0)] for phi, ips, mini in zip(phi_lep, sip3d_lep, mini_lep)])
                dphi_lep_met = np.array([np.array([phi - m_phi for phi in leps])
                                            for leps, m_phi in zip(phi_mini_lep, met_phi)])

                dphi_lep_met = np.array([ np.array([phi + 2*np.pi if phi < -np.pi else phi for phi in leps]) for leps in dphi_lep_met])
                dphi_lep_met = np.array([ np.array([phi - 2*np.pi if phi >= np.pi else phi for phi in leps]) for leps in dphi_lep_met])
                dphi_lep_met = np.abs(dphi_lep_met) 
                print '\ncreating masks and weights'
                print '-> bjet masks'

                print '-> lepton masks'
                only_2_golden_lep = np.array([True if len(lep) == 2 and two else False for lep, two in zip(pt_35_lep, two_lep)])
                only_1_golden_lep = np.array([True if len(lep) == 1 and one else False for lep, one in zip(pt_35_lep, one_lep)])

                only_1_golden_less12_lep = np.array([True if len(lep[lep<12]) == 1 and one else False for lep, one in zip(pt_35_lep, one_lep)])
                only_1_golden_less10_lep = np.array([True if len(lep[lep<10]) == 1 and one else False for lep, one in zip(pt_35_lep, one_lep)])
                only_1_golden_less7_lep = np.array([True if len(lep[lep<7]) == 1 and one else False for lep, one in zip(pt_35_lep, one_lep)])

                only_1_goldplus_lep = np.array([True if len(lep) == 1 and one and len(charge[charge>0])>0 else False for lep, one, charge in zip(pt_35_lep, one_lep, ch_35_lep)])
                only_1_golden_el = np.array([True if len(lep) == 1 and one else False for lep, one in zip(pt_5_el, one_lep)])
                only_1_golden_mu = np.array([True if len(lep) == 1 and one else False for lep, one in zip(pt_35_mu, one_lep)])
#                only_1_less_pt12_lep = np.array([True if len(lep) == 1 and not two_leps else False for lep, two_leps in zip(pt_less_12_lep, only_2_lep)])

                only_2_bronze_lep = np.array([True if two and not gold else False for gold, two in zip(only_2_golden_lep, two_lep)])
                only_1_bronze_lep = np.array([True if one and not gold else False for gold, one in zip(only_1_golden_lep, one_lep)])
                only_1_bronze_el = np.array([True if one and not gold else False for gold, one in zip(only_1_golden_el, one_lep)])
                only_1_bronze_mu = np.array([True if one and not gold else False for gold, one in zip(only_1_golden_mu, one_lep)])

                dphi_2lep_met = np.array([ dphi if lep else np.array([None, None]) for dphi, lep in zip(dphi_lep_met, only_2_golden_lep)])
                dphi_lep_met = np.array([ dphi if lep else np.array([None]) for dphi, lep in zip(dphi_lep_met, only_1_golden_lep)])

                dphi_g1p5 = np.array([True if np.all(lep>1.5) else False for lep in dphi_lep_met])
                dphi_l1p5 = np.array([True if np.all(lep<1.5) else False for lep in dphi_lep_met])
                dphi_2_g1p5 = np.array([True if np.all(lep>1.5) else False for lep in dphi_2lep_met])
                dphi_2_l1p5 = np.array([True if np.all(lep<1.5) else False for lep in dphi_2lep_met])

                dphi_lep_s_l0p3 = np.array([True if np.all(lep<0.3) else False for lep in dphi_lep_met])

                met_200 = met > 200

                risr_0p95 = risr > 0.95
                risr_0p9 = risr > 0.9
                risr_0p8 = risr > 0.8
                risr_lp9 = risr < 0.9
                risr_lp95 = risr < 0.95

                ptisr_200 = ptisr > 200

                has_mscut1 = ms < 80.
                has_mscut2a = ms > 80.
                has_mscut2b = ms < 160.
                has_mscut2 = np.all([has_mscut2a, has_mscut2b], axis=0)
                has_mscut3 = ms > 160.

                has_dphilep = dphi_lep_s_l0p3 > 0.

#                for sv, eta in zip(nsv_eta_s, eta_sv_s):
#                    print sv, eta
                print 'incrementing tables'
                has_greateq_1sv = nsv_eta >= 1
                has_ge1_s_sv = nsv_eta_s >= 1
                has_ge1_isrv = nsv_eta_isr >= 1
                has_0st = njet_s + nsv_eta_s < 1
                has_greateq_1st = njet_s + nsv_eta_s >= 1
                has_greateq_2st = njet_s + nsv_eta_s >= 2
                has_greateq_2sv = nsv_eta_s >= 2
                has_greateq_2isrv = nsv_eta_isr >= 2
                has_greateq_1sb = nbjet_s >= 1
                e_ge1sv = np.all([met_200, risr_0p8, ptisr_200, has_greateq_1sv], axis=0)
                e_ge2sv = np.all([met_200, risr_0p8, ptisr_200, has_greateq_2sv], axis=0)
                e_ge1sb = np.all([met_200, risr_0p8, ptisr_200, has_greateq_1sb], axis=0)
                two_l = np.all([met_200, risr_0p8, ptisr_200, only_2_golden_lep], axis=0)
                two_l_ge1st = np.all([met_200, risr_0p8, ptisr_200, only_2_golden_lep, has_greateq_1st], axis=0)
                two_l_ge1sv = np.all([met_200, risr_0p8, ptisr_200, only_2_golden_lep, has_greateq_1sv], axis=0)
                one_l = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep], axis=0)
                one_el = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_el], axis=0)
                one_mu = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_mu], axis=0)
                has_no_medium = nbjet_s < 1
#                one_0sb = np.all([met_200, risr_0p8, ptisr_200, only_1_less_pt12_lep, has_no_medium], axis=0)
                one_0sb = np.all([met_200, risr_0p8, ptisr_200, has_no_medium], axis=0)
                one_l_greateq_1st = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_greateq_1st], axis=0)
                one_el_greateq_1st = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_el, has_greateq_1st], axis=0)
                one_mu_greateq_1st = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_mu, has_greateq_1st], axis=0)
                has_1_medium = nbjet_s == 1
                one_1sb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_1_medium], axis=0)
#                has_greateq_1sv = nsv_eta_s >= 1
                one_ge1v = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_greateq_1sv], axis=0)
                has_0_s_jets = njet_s < 1
                has_greateq_1_s_jets = njet_s >= 1
                one_ge1sb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_greateq_1sb], axis=0)
                one_ge1sj = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_greateq_1_s_jets], axis=0)
                has_0_s_sv = nsv_eta_s < 1
                one_l_0ssv = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_s_sv], axis=0)
                has_0_isr_medium = nbjet_isr < 1
                has_0_isr_t = nbjet_isr + nsv_eta_isr < 1
                has_ge1_isr_medium = nbjet_isr >= 1
                has_ge2_isr_medium = nbjet_isr >= 2
                has_1_isr_medium = nbjet_isr == 1
                has_ge1_s_medium = nbjet_s >= 1
                has_greateq_1_isr_t = nbjet_isr + nsv_eta_isr >= 1
                e_isr_ge1v = np.all([met_200, risr_0p8, ptisr_200, has_ge1_isrv], axis=0)
                e_isr_ge2v = np.all([met_200, risr_0p8, ptisr_200, has_greateq_2isrv], axis=0)
                e_isr_ge1b = np.all([met_200, risr_0p8, ptisr_200, has_ge1_isr_medium], axis=0)
                two_l_0isrb = np.all([met_200, risr_0p8, ptisr_200, only_2_golden_lep, has_0_isr_medium], axis=0)
                two_l_0isrb_ge1st = np.all([met_200, risr_0p8, ptisr_200, only_2_golden_lep, has_0_isr_medium, has_greateq_1st], axis=0)
                two_l_0isrb_ge1sv = np.all([met_200, risr_0p8, ptisr_200, only_2_golden_lep, has_0_isr_medium, has_greateq_1sv], axis=0)
                two_l_0isrb_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, only_2_golden_lep, has_0_isr_medium], axis=0)
                two_l_0isrb_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, only_2_golden_lep, has_0_isr_medium], axis=0)
                two_l_0isrb_r95 = np.all([met_200, risr_0p95, ptisr_200, only_2_golden_lep, has_0_isr_medium], axis=0)
                two_l_0isrb_ge1st_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, only_2_golden_lep, has_0_isr_medium, has_greateq_1st], axis=0)
                two_l_0isrb_ge1st_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, only_2_golden_lep, has_0_isr_medium, has_greateq_1st], axis=0)
                two_l_0isrb_ge1st_r95 = np.all([met_200, risr_0p95, ptisr_200, only_2_golden_lep, has_0_isr_medium, has_greateq_1st], axis=0)
                two_l_0isrb_ge1sv_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, only_2_golden_lep, has_0_isr_medium, has_greateq_1sv], axis=0)
                two_l_0isrb_ge1sv_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, only_2_golden_lep, has_0_isr_medium, has_greateq_1sv], axis=0)
                two_l_0isrb_ge1sv_r95 = np.all([met_200, risr_0p95, ptisr_200, only_2_golden_lep, has_0_isr_medium, has_greateq_1sv], axis=0)
                two_l_ge1isrb = np.all([met_200, risr_0p8, ptisr_200, only_2_golden_lep, has_ge1_isr_medium], axis=0)
                two_l_ge1isrb_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, only_2_golden_lep, has_ge1_isr_medium], axis=0)
                two_l_ge1isrb_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, only_2_golden_lep, has_ge1_isr_medium], axis=0)
                two_l_ge1isrb_r95 = np.all([met_200, risr_0p95, ptisr_200, only_2_golden_lep, has_ge1_isr_medium], axis=0)
                one_l_0isrb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium], axis=0)
                one_l_0isrt = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_t], axis=0)
                one_l_0isrb_0ssv = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_0_s_sv], axis=0)
                one_l_0isrb_0s_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, only_1_golden_lep, has_0st], axis=0)
                one_l_0isrb_0s_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, only_1_golden_lep, has_0st], axis=0)
                one_l_0isrb_0s_r95 = np.all([met_200, risr_0p95, ptisr_200, only_1_golden_lep, has_0st], axis=0)
                one_l_1isrb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_1_isr_medium], axis=0)
                one_l_greateq_1isrb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_ge1_isr_medium], axis=0)
                one_l_greateq_1isrt = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_greateq_1_isr_t], axis=0)
                one_l_ge2isrb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_ge2_isr_medium], axis=0)
                one_l_ge2isrb_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, only_1_golden_lep, has_ge2_isr_medium], axis=0)
                one_l_ge2isrb_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, only_1_golden_lep, has_ge2_isr_medium], axis=0)
                one_l_ge2isrb_r95 = np.all([met_200, risr_0p95, ptisr_200, only_1_golden_lep, has_ge2_isr_medium], axis=0)
                one_el_0isrb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_el, has_0_isr_medium], axis=0)
                one_mu_0isrb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_mu, has_0_isr_medium], axis=0)
#                one_l_greateq_1isrb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_ge1_isr_medium], axis=0)
#                one_l_greateq_1isrt = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_greateq_1_isr_t], axis=0)
#                one_l_1isrb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_1_isr_medium], axis=0)
                one_l_0isrb_ge1sj_m1_r95 = np.all([met_200, risr_0p95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1_s_jets, has_mscut1], axis=0)
#                one_l_0isrb_ge1sj_mscut2 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1_s_jets, has_mscut2], axis=0)
                one_l_0isrb_ge1sj = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1_s_jets], axis=0)
                one_l_0isrb_ge1sb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_medium], axis=0)
                one_l_0isrb_ge1sv = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv], axis=0)
                one_l_0isrb_ge1v = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1sv], axis=0)
                one_l_0isrb_ge1sb_ge1sv = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_ge1_s_medium], axis=0)
                one_l_0isrb_ge1sb_ge1sv_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_ge1_s_medium], axis=0)
                one_l_0isrb_ge1sb_ge1sv_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_ge1_s_medium], axis=0)
                one_l_0isrb_ge1sb_ge1sv_r95 = np.all([met_200, risr_0p95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_ge1_s_medium], axis=0)
#############
#                one_l_0isrb_ge1sv_philmetg = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv], axis=0)
                one_l_0isrb_ge1sv_philmetg = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, dphi_g1p5], axis=0)
                one_l_0isrb_ge1sv_philmetl = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, dphi_l1p5], axis=0)
#                one_l_0isrb_ge1sv_philmetl = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv], axis=0)
                one_l_0isrb_ge1st_m1_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st, has_mscut1], axis=0)
                one_l_0isrb_ge1st_m1_r95 = np.all([met_200, risr_0p95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st, has_mscut1], axis=0)
                one_l_0isrb_ge1st_mscut2 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st, has_mscut2], axis=0)
                one_l_0isrb_ge1st_mscut3 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st, has_mscut3], axis=0)
                one_l_0isrb_ge1st_m3_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st, has_mscut3], axis=0)
                one_l_0isrb_ge1v_mscut1 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1sv, has_mscut1], axis=0)
                one_l_0isrb_ge1sv_mscut1 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut1, dphi_l1p5], axis=0)
                one_l_0isrb_ge1sv_m1_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut1], axis=0)
                one_l_0isrb_ge1sv_m1_r95 = np.all([met_200, risr_0p95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut1], axis=0)
                one_l_0isrb_ge1sv_mscut1_philmetl = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut1, dphi_l1p5], axis=0)
                one_l_0isrb_ge1sv_mscut3_philmetg = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut3, dphi_g1p5], axis=0)
                one_l_0isrb_ge1sv_mscut1_philmetl_r95 = np.all([met_200, risr_0p95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut1, dphi_l1p5], axis=0)
                one_l_0isrb_ge1sv_m3_philmetg_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut3, dphi_g1p5], axis=0)
#                one_l_0isrb_ge1sv_mscut1_philmetl = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut1], axis=0)
#                one_l_0isrb_ge1sv_mscut3_philmetg = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut3], axis=0)
                one_l_0isrb_ge2sv = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_2sv], axis=0)
                one_l_0isrb_ge2st = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_2st], axis=0)
                one_l_0isrb_ge2st_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_2st], axis=0)
                one_l_0isrb_ge2st_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_2st], axis=0)
                one_l_0isrb_ge2st_r95 = np.all([met_200, risr_0p95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_2st], axis=0)
#                one_l_0isrb_ge1sj_mscut3 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1_s_jets, has_mscut3], axis=0)
                one_l_0isrb_ge1sj_mscut1 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1_s_jets, has_mscut1], axis=0)
                one_l_0isrb_ge1sj_mscut2 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1_s_jets, has_mscut2], axis=0)
                one_l_0isrb_ge1sb_mscut1 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_medium, has_mscut1], axis=0)
                one_l_0isrb_ge1sb_m1_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_medium, has_mscut1], axis=0)
                one_l_0isrb_ge1sb_m1_r95 = np.all([met_200, risr_0p95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_medium, has_mscut1], axis=0)
                one_l_0isrb_ge1sb_mscut2 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_medium, has_mscut2], axis=0)
                one_l_0isrb_ge1sb_mscut3 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_medium, has_mscut3], axis=0)
                one_l_0isrb_ge1sb_m3_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_medium, has_mscut3], axis=0)
                one_l_0isrb_ge1st = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st], axis=0)
                one_l_0isrb_ge1st_mscut1 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st, has_mscut1], axis=0)
                one_l_0isrb_ge1st_mscut2 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st, has_mscut2], axis=0)
                one_l_0isrb_ge1st_mscut3 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st, has_mscut3], axis=0)
                one_l_0isrb_ge1sv_mscut1 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut1], axis=0)
                one_l_0isrb_ge1sv_mscut1_philmetl = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut1, dphi_l1p5], axis=0)
                one_l_0isrb_ge1sv_m1_philmetl_r95 = np.all([met_200, risr_0p95, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut1, dphi_l1p5], axis=0)
                one_l_0isrb_ge1sv_mscut1_philmetg = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut1, dphi_g1p5], axis=0)
                one_l_0isrb_ge1sv_mscut2 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut2], axis=0)
                one_l_0isrb_ge1sv_mscut3 = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut3], axis=0)
                one_l_0isrb_ge1sv_mscut3_philmetg = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut3, dphi_g1p5], axis=0)
                one_l_0isrb_ge1sv_mscut3_philmetg = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_ge1_s_sv, has_mscut3, dphi_g1p5], axis=0)
                one_l_0isrb_ge1st_philmetg = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st, dphi_g1p5], axis=0)
                one_l_0isrb_ge1st_philmetl = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st, dphi_l1p5], axis=0)
                one_l_0isrb_ge1sv_pt12_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_less12_lep], axis=0)
                one_l_0isrb_ge1sv_pt12_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_less12_lep], axis=0)
                one_l_0isrb_ge1sv_pt12_r95 = np.all([met_200, risr_0p95, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_less12_lep], axis=0)
                one_l_0isrb_ge1sv_pt10_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_less10_lep], axis=0)
                one_l_0isrb_ge1sv_pt10_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_less10_lep], axis=0)
                one_l_0isrb_ge1sv_pt10_r95 = np.all([met_200, risr_0p95, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_less10_lep], axis=0)
                one_l_0isrb_ge1sv_pt7_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_less7_lep], axis=0)
                one_l_0isrb_ge1sv_pt7_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_less7_lep], axis=0)
                one_l_0isrb_ge1sv_pt7_r95 = np.all([met_200, risr_0p95, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_less7_lep], axis=0)
                one_l_0isrb_ge1sv_dphlep_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_lep, has_dphilep], axis=0)
                one_l_0isrb_ge1sv_dphlep_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_lep, has_dphilep], axis=0)
                one_l_0isrb_ge1sv_dphlep_r95 = np.all([met_200, risr_0p95, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_lep, has_dphilep], axis=0)
                one_l_0isrb_ge1sv_dphlep_pt7_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_less7_lep, has_dphilep], axis=0)
                one_l_0isrb_ge1sv_dphlep_pt7_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_less7_lep, has_dphilep], axis=0)
                one_l_0isrb_ge1sv_dphlep_pt7_r95 = np.all([met_200, risr_0p95, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_golden_less7_lep, has_dphilep], axis=0)
                one_l_0isrb_ge1sv_chp_r8 = np.all([met_200, risr_0p8, risr_lp9, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_goldplus_lep], axis=0)
                one_l_0isrb_ge1sv_chp_r9 = np.all([met_200, risr_0p9, risr_lp95, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_goldplus_lep], axis=0)
                one_l_0isrb_ge1sv_chp_r95 = np.all([met_200, risr_0p95, ptisr_200, has_0_isr_medium, has_ge1_s_sv, only_1_goldplus_lep], axis=0)
#                one_l_0isrb_ge1st_philmetg = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st], axis=0)
#                one_l_0isrb_ge1st_philmetl = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_lep, has_0_isr_medium, has_greateq_1st], axis=0)
                one_el_0isrb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_el, has_0_isr_medium], axis=0)
                one_el_0isrb_ge1st = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_el, has_0_isr_medium, has_greateq_1st], axis=0)
                one_el_0isrb_ge1sj = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_el, has_0_isr_medium, has_greateq_1_s_jets], axis=0)
                one_mu_0isrb = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_mu, has_0_isr_medium], axis=0)
                one_mu_0isrb_ge1st = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_mu, has_0_isr_medium, has_greateq_1st], axis=0)
                one_mu_0isrb_ge1sj = np.all([met_200, risr_0p8, ptisr_200, only_1_golden_mu, has_0_isr_medium, has_greateq_1_s_jets], axis=0)

                # s-style regions
                hist[sample][tree_name][reference.index('total')] += len(met)
                hist[sample][tree_name][reference.index('ge1sv')] += len(met[e_ge1sv])
                hist[sample][tree_name][reference.index('ge2sv')] += len(met[e_ge2sv])
                hist[sample][tree_name][reference.index('ge1sb')] += len(met[e_ge1sb])
                hist[sample][tree_name][reference.index('2l')] += len(met[two_l])
                hist[sample][tree_name][reference.index('2l_ge1st')] += len(met[two_l_ge1st])
                hist[sample][tree_name][reference.index('2l_ge1sv')] += len(met[two_l_ge1sv])
                hist[sample][tree_name][reference.index('1l')] += len(met[one_l])
                hist[sample][tree_name][reference.index('1_el')] += len(met[one_el])
                hist[sample][tree_name][reference.index('1_mu')] += len(met[one_mu])
                hist[sample][tree_name][reference.index('1l_0sb')] += len(met[one_0sb])
                hist[sample][tree_name][reference.index('1l_ge1st')] += len(met[one_l_greateq_1st])
                hist[sample][tree_name][reference.index('1_el_ge1st')] += len(met[one_el_greateq_1st])
                hist[sample][tree_name][reference.index('1_mu_ge1st')] += len(met[one_mu_greateq_1st])
                hist[sample][tree_name][reference.index('1l_1sb')] += len(met[one_1sb])
                hist[sample][tree_name][reference.index('1l_ge1v')] += len(met[one_ge1v])
                hist[sample][tree_name][reference.index('1l_ge1sb')] += len(met[one_ge1sb])
                hist[sample][tree_name][reference.index('1l_ge1sj')] += len(met[one_ge1sj])
                hist[sample][tree_name][reference.index('1l_0ssv')] += len(met[one_l_0ssv])
                
                # isr-style regions
                hist[sample][tree_name][reference.index('total')] += len(met)
                hist[sample][tree_name][reference.index('isr_ge1v')] += len(met[e_isr_ge1v])
                hist[sample][tree_name][reference.index('isr_ge2v')] += len(met[e_isr_ge2v])
                hist[sample][tree_name][reference.index('isr_ge1b')] += len(met[e_isr_ge1b])
                hist[sample][tree_name][reference.index('2l_0isrb')] += len(met[two_l_0isrb])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1st')] += len(met[two_l_0isrb_ge1st])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1sv')] += len(met[two_l_0isrb_ge1sv])
                hist[sample][tree_name][reference.index('2l_ge1isrb')] += len(met[two_l_ge1isrb])
                hist[sample][tree_name][reference.index('2l_0isrb_r8')] += len(met[two_l_0isrb_r8])
                hist[sample][tree_name][reference.index('2l_0isrb_r9')] += len(met[two_l_0isrb_r9])
                hist[sample][tree_name][reference.index('2l_0isrb_r95')] += len(met[two_l_0isrb_r95])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1st_r8')] += len(met[two_l_0isrb_ge1st_r8])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1st_r9')] += len(met[two_l_0isrb_ge1st_r9])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1st_r95')] += len(met[two_l_0isrb_ge1st_r95])
                hist[sample][tree_name][reference.index('2l_ge1isrb_r8')] += len(met[two_l_ge1isrb_r8])
                hist[sample][tree_name][reference.index('2l_ge1isrb_r9')] += len(met[two_l_ge1isrb_r9])
                hist[sample][tree_name][reference.index('2l_ge1isrb_r95')] += len(met[two_l_ge1isrb_r95])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1sv_r8')] += len(met[two_l_0isrb_ge1sv_r8])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1sv_r9')] += len(met[two_l_0isrb_ge1sv_r9])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1sv_r95')] += len(met[two_l_0isrb_ge1sv_r95])

                hist[sample][tree_name][reference.index('1l_0isrb')] += len(met[one_l_0isrb])
                hist[sample][tree_name][reference.index('1l_0isrb_0ssv')] += len(met[one_l_0isrb_0ssv])
                hist[sample][tree_name][reference.index('1l_0isrb_0s_r8')] += len(met[one_l_0isrb_0s_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_0s_r9')] += len(met[one_l_0isrb_0s_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_0s_r95')] += len(met[one_l_0isrb_0s_r95])
                hist[sample][tree_name][reference.index('1l_ge2isrb')] += len(met[one_l_ge2isrb])
                hist[sample][tree_name][reference.index('1l_ge2isrb_r8')] += len(met[one_l_ge2isrb_r8])
                hist[sample][tree_name][reference.index('1l_ge2isrb_r9')] += len(met[one_l_ge2isrb_r9])
                hist[sample][tree_name][reference.index('1l_ge2isrb_r95')] += len(met[one_l_ge2isrb_r95])
                hist[sample][tree_name][reference.index('1_el_0isrb')] += len(met[one_el_0isrb])
                hist[sample][tree_name][reference.index('1_mu_0isrb')] += len(met[one_mu_0isrb])
                hist[sample][tree_name][reference.index('1l_ge1isrb')] += len(met[one_l_greateq_1isrb])
                hist[sample][tree_name][reference.index('1l_ge1isrt')] += len(met[one_l_greateq_1isrt])
                hist[sample][tree_name][reference.index('1l_1isrb')] += len(met[one_l_1isrb])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sj_mscut1')] += len(met[one_l_0isrb_ge1sj_mscut1])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sj_m1_r95')] += len(met[one_l_0isrb_ge1sj_m1_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sj_mscut2')] += len(met[one_l_0isrb_ge1sj_mscut2])
                hist[sample][tree_name][reference.index('1_el_0isrb_ge1sj')] += len(met[one_el_0isrb_ge1sj])
                hist[sample][tree_name][reference.index('1_mu_0isrb_ge1sj')] += len(met[one_mu_0isrb_ge1sj])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sj')] += len(met[one_l_0isrb_ge1sj])
                hist[sample][tree_name][reference.index('1l_0isrt')] += len(met[one_l_0isrt])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb')] += len(met[one_l_0isrb_ge1sb])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv')] += len(met[one_l_0isrb_ge1sv])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1v')] += len(met[one_l_0isrb_ge1v])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_ge1sv')] += len(met[one_l_0isrb_ge1sb_ge1sv])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_ge1sv_r8')] += len(met[one_l_0isrb_ge1sb_ge1sv_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_ge1sv_r9')] += len(met[one_l_0isrb_ge1sb_ge1sv_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_ge1sv_r95')] += len(met[one_l_0isrb_ge1sb_ge1sv_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_philmetg')] += len(met[one_l_0isrb_ge1sv_philmetg])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_philmetl')] += len(met[one_l_0isrb_ge1sv_philmetl])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st')] += len(met[one_l_0isrb_ge1st])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_mscut1')] += len(met[one_l_0isrb_ge1st_mscut1])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_m1_r9')] += len(met[one_l_0isrb_ge1st_m1_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_m1_r95')] += len(met[one_l_0isrb_ge1st_m1_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_mscut2')] += len(met[one_l_0isrb_ge1st_mscut2])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_mscut3')] += len(met[one_l_0isrb_ge1st_mscut3])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_m3_r8')] += len(met[one_l_0isrb_ge1st_m3_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_mscut1')] += len(met[one_l_0isrb_ge1sv_mscut1])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_m1_r95')] += len(met[one_l_0isrb_ge1sv_m1_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_m1_r9')] += len(met[one_l_0isrb_ge1sv_m1_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_mscut2')] += len(met[one_l_0isrb_ge1sv_mscut2])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_mscut3')] += len(met[one_l_0isrb_ge1sv_mscut3])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1v_mscut1')] += len(met[one_l_0isrb_ge1v_mscut1])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_mscut1_philmetl')] += len(met[one_l_0isrb_ge1sv_mscut1_philmetl])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_m1_philmetl_r95')] += len(met[one_l_0isrb_ge1sv_m1_philmetl_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_mscut1_philmetg')] += len(met[one_l_0isrb_ge1sv_mscut1_philmetg])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_mscut3_philmetg')] += len(met[one_l_0isrb_ge1sv_mscut3_philmetg])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_m3_philmetg_r8')] += len(met[one_l_0isrb_ge1sv_m3_philmetg_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge2sv')] += len(met[one_l_0isrb_ge2sv])
                hist[sample][tree_name][reference.index('1l_0isrb_ge2st_r9')] += len(met[one_l_0isrb_ge2st_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge2st_r8')] += len(met[one_l_0isrb_ge2st_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge2st_r95')] += len(met[one_l_0isrb_ge2st_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_m1_r9')] += len(met[one_l_0isrb_ge1sb_m1_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_m1_r95')] += len(met[one_l_0isrb_ge1sb_m1_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_mscut2')] += len(met[one_l_0isrb_ge1sb_mscut2])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_m3_r8')] += len(met[one_l_0isrb_ge1sb_m3_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st')] += len(met[one_l_0isrb_ge1st])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_philmetg')] += len(met[one_l_0isrb_ge1st_philmetg])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_philmetl')] += len(met[one_l_0isrb_ge1st_philmetl])
                hist[sample][tree_name][reference.index('1_el_0isrb_ge1st')] += len(met[one_el_0isrb_ge1st])
                hist[sample][tree_name][reference.index('1_mu_0isrb_ge1st')] += len(met[one_mu_0isrb_ge1st])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt12_r8')] += len(met[one_l_0isrb_ge1sv_pt12_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt12_r9')] += len(met[one_l_0isrb_ge1sv_pt12_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt12_r95')] += len(met[one_l_0isrb_ge1sv_pt12_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt10_r8')] += len(met[one_l_0isrb_ge1sv_pt10_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt10_r9')] += len(met[one_l_0isrb_ge1sv_pt10_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt10_r95')] += len(met[one_l_0isrb_ge1sv_pt10_r95])                
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt7_r8')] += len(met[one_l_0isrb_ge1sv_pt7_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt7_r9')] += len(met[one_l_0isrb_ge1sv_pt7_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt7_r95')] += len(met[one_l_0isrb_ge1sv_pt7_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_dphlep_r8')] += len(met[one_l_0isrb_ge1sv_dphlep_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_dphlep_r9')] += len(met[one_l_0isrb_ge1sv_dphlep_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_dphlep_r95')] += len(met[one_l_0isrb_ge1sv_dphlep_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_dphlep_pt7_r8')] += len(met[one_l_0isrb_ge1sv_dphlep_pt7_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_dphlep_pt7_r9')] += len(met[one_l_0isrb_ge1sv_dphlep_pt7_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_dphlep_pt7_r95')] += len(met[one_l_0isrb_ge1sv_dphlep_pt7_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_chp_r8')] += len(met[one_l_0isrb_ge1sv_chp_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_chp_r9')] += len(met[one_l_0isrb_ge1sv_chp_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_chp_r95')] += len(met[one_l_0isrb_ge1sv_chp_r95])

                # weighted entries
                # s-style regions
                hist[sample][tree_name][reference.index('total_w')] += np.sum(weight) 
                hist[sample][tree_name][reference.index('ge1sv_w')] += np.sum(weight[e_ge1sv])
                hist[sample][tree_name][reference.index('ge2sv_w')] += np.sum(weight[e_ge2sv])
                hist[sample][tree_name][reference.index('ge1sb_w')] += np.sum(weight[e_ge1sb])
                hist[sample][tree_name][reference.index('2l_w')] += np.sum(weight[two_l])
                hist[sample][tree_name][reference.index('2l_ge1st_w')] += np.sum(weight[two_l_ge1st])
                hist[sample][tree_name][reference.index('2l_ge1sv_w')] += np.sum(weight[two_l_ge1sv])
                hist[sample][tree_name][reference.index('2l_0isrb_r8_w')] += np.sum(weight[two_l_0isrb_r8])
                hist[sample][tree_name][reference.index('2l_0isrb_r9_w')] += np.sum(weight[two_l_0isrb_r9])
                hist[sample][tree_name][reference.index('2l_0isrb_r95_w')] += np.sum(weight[two_l_0isrb_r95])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1st_r8_w')] += np.sum(weight[two_l_0isrb_ge1st_r8])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1st_r9_w')] += np.sum(weight[two_l_0isrb_ge1st_r9])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1st_r95_w')] += np.sum(weight[two_l_0isrb_ge1st_r95])
                hist[sample][tree_name][reference.index('2l_ge1isrb_r8_w')] += np.sum(weight[two_l_ge1isrb_r8])
                hist[sample][tree_name][reference.index('2l_ge1isrb_r9_w')] += np.sum(weight[two_l_ge1isrb_r9])
                hist[sample][tree_name][reference.index('2l_ge1isrb_r95_w')] += np.sum(weight[two_l_ge1isrb_r95])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1sv_r8_w')] += np.sum(weight[two_l_0isrb_ge1sv_r8])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1sv_r9_w')] += np.sum(weight[two_l_0isrb_ge1sv_r9])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1sv_r95_w')] += np.sum(weight[two_l_0isrb_ge1sv_r95])
                hist[sample][tree_name][reference.index('2l_0isrb_w')] += np.sum(weight[two_l_0isrb])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1st_w')] += np.sum(weight[two_l_0isrb_ge1st])
                hist[sample][tree_name][reference.index('2l_0isrb_ge1sv_w')] += np.sum(weight[two_l_0isrb_ge1sv])
                hist[sample][tree_name][reference.index('2l_ge1isrb_w')] += np.sum(weight[two_l_ge1isrb])


                hist[sample][tree_name][reference.index('1l_w')] += np.sum(weight[one_l])
                hist[sample][tree_name][reference.index('1_el_w')] += np.sum(weight[one_el])
                hist[sample][tree_name][reference.index('1_mu_w')] += np.sum(weight[one_mu])
                hist[sample][tree_name][reference.index('1l_0sb_w')] += np.sum(weight[one_0sb])
                hist[sample][tree_name][reference.index('1l_ge1st_w')] += np.sum(weight[one_l_greateq_1st])
                hist[sample][tree_name][reference.index('1_el_ge1st_w')] += np.sum(weight[one_el_greateq_1st])
                hist[sample][tree_name][reference.index('1_mu_ge1st_w')] += np.sum(weight[one_mu_greateq_1st])
                hist[sample][tree_name][reference.index('1l_1sb_w')] += np.sum(weight[one_1sb])
                hist[sample][tree_name][reference.index('1l_ge1v_w')] += np.sum(weight[one_ge1v])
                hist[sample][tree_name][reference.index('1l_ge1sb_w')] += np.sum(weight[one_ge1sb])
                hist[sample][tree_name][reference.index('1l_ge1sj_w')] += np.sum(weight[one_ge1sj])
                hist[sample][tree_name][reference.index('1l_0ssv_w')] += np.sum(weight[one_l_0ssv])
               
                hist[sample][tree_name][reference.index('total_w')] += np.sum(weight) 
                hist[sample][tree_name][reference.index('isr_ge1v_w')] += np.sum(weight[e_isr_ge1v])
                hist[sample][tree_name][reference.index('isr_ge2v_w')] += np.sum(weight[e_isr_ge2v])
                hist[sample][tree_name][reference.index('isr_ge1b_w')] += np.sum(weight[e_isr_ge1b])
                hist[sample][tree_name][reference.index('1l_0isrb_w')] += np.sum(weight[one_l_0isrb])
                hist[sample][tree_name][reference.index('1l_0isrb_0ssv_w')] += np.sum(weight[one_l_0isrb_0ssv])
                hist[sample][tree_name][reference.index('1l_0isrb_0s_r8_w')] += np.sum(weight[one_l_0isrb_0s_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_0s_r9_w')] += np.sum(weight[one_l_0isrb_0s_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_0s_r95_w')] += np.sum(weight[one_l_0isrb_0s_r95])
                hist[sample][tree_name][reference.index('1l_ge2isrb_w')] += np.sum(weight[one_l_ge2isrb])
                hist[sample][tree_name][reference.index('1l_ge2isrb_r8_w')] += np.sum(weight[one_l_ge2isrb_r8])
                hist[sample][tree_name][reference.index('1l_ge2isrb_r9_w')] += np.sum(weight[one_l_ge2isrb_r9])
                hist[sample][tree_name][reference.index('1l_ge2isrb_r95_w')] += np.sum(weight[one_l_ge2isrb_r95])
                hist[sample][tree_name][reference.index('1_el_0isrb_w')] += np.sum(weight[one_el_0isrb])
                hist[sample][tree_name][reference.index('1_mu_0isrb_w')] += np.sum(weight[one_mu_0isrb])
                hist[sample][tree_name][reference.index('1l_ge1isrb_w')] += np.sum(weight[one_l_greateq_1isrb])
                hist[sample][tree_name][reference.index('1l_ge1isrt_w')] += np.sum(weight[one_l_greateq_1isrt])
                hist[sample][tree_name][reference.index('1l_1isrb_w')] += np.sum(weight[one_l_1isrb])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sj_m1_r95_w')] += np.sum(weight[one_l_0isrb_ge1sj_m1_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sj_mscut1_w')] += np.sum(weight[one_l_0isrb_ge1sj_mscut1])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sj_mscut2_w')] += np.sum(weight[one_l_0isrb_ge1sj_mscut2])
                hist[sample][tree_name][reference.index('1_el_0isrb_ge1sj_w')] += np.sum(weight[one_el_0isrb_ge1sj])
                hist[sample][tree_name][reference.index('1_mu_0isrb_ge1sj_w')] += np.sum(weight[one_mu_0isrb_ge1sj])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sj_w')] += np.sum(weight[one_l_0isrb_ge1sj])
                hist[sample][tree_name][reference.index('1l_0isrt_w')] += np.sum(weight[one_l_0isrt])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_w')] += np.sum(weight[one_l_0isrb_ge1sb])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_w')] += np.sum(weight[one_l_0isrb_ge1sv])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1v_w')] += np.sum(weight[one_l_0isrb_ge1v])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_ge1sv_w')] += np.sum(weight[one_l_0isrb_ge1sb_ge1sv])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_ge1sv_r8_w')] += np.sum(weight[one_l_0isrb_ge1sb_ge1sv_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_ge1sv_r9_w')] += np.sum(weight[one_l_0isrb_ge1sb_ge1sv_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_ge1sv_r95_w')] += np.sum(weight[one_l_0isrb_ge1sb_ge1sv_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_philmetg_w')] += np.sum(weight[one_l_0isrb_ge1sv_philmetg])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_philmetl_w')] += np.sum(weight[one_l_0isrb_ge1sv_philmetl])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_mscut1_w')] += np.sum(weight[one_l_0isrb_ge1st_mscut1])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_m3_r8_w')] += np.sum(weight[one_l_0isrb_ge1st_m3_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_m1_r9_w')] += np.sum(weight[one_l_0isrb_ge1st_m1_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_m1_r95_w')] += np.sum(weight[one_l_0isrb_ge1st_m1_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_mscut2_w')] += np.sum(weight[one_l_0isrb_ge1st_mscut2])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_mscut3_w')] += np.sum(weight[one_l_0isrb_ge1st_mscut3])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_mscut1_w')] += np.sum(weight[one_l_0isrb_ge1sv_mscut1])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_m1_r9_w')] += np.sum(weight[one_l_0isrb_ge1sv_m1_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_m1_r95_w')] += np.sum(weight[one_l_0isrb_ge1sv_m1_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1v_mscut1_w')] += np.sum(weight[one_l_0isrb_ge1v_mscut1])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_mscut3_w')] += np.sum(weight[one_l_0isrb_ge1sv_mscut3])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_m1_philmetl_r95_w')] += np.sum(weight[one_l_0isrb_ge1sv_m1_philmetl_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_m3_philmetg_r8_w')] += np.sum(weight[one_l_0isrb_ge1sv_m3_philmetg_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge2sv_w')] += np.sum(weight[one_l_0isrb_ge2sv])
                hist[sample][tree_name][reference.index('1l_0isrb_ge2st_w')] += np.sum(weight[one_l_0isrb_ge2st])
                hist[sample][tree_name][reference.index('1l_0isrb_ge2st_r8_w')] += np.sum(weight[one_l_0isrb_ge2st_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge2st_r9_w')] += np.sum(weight[one_l_0isrb_ge2st_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge2st_r95_w')] += np.sum(weight[one_l_0isrb_ge2st_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_m1_r9_w')] += np.sum(weight[one_l_0isrb_ge1sb_m1_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_m1_r95_w')] += np.sum(weight[one_l_0isrb_ge1sb_m1_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_mscut2_w')] += np.sum(weight[one_l_0isrb_ge1sb_mscut2])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sb_m3_r8_w')] += np.sum(weight[one_l_0isrb_ge1sb_m3_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_w')] += np.sum(weight[one_l_0isrb_ge1st])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_philmetg_w')] += np.sum(weight[one_l_0isrb_ge1st_philmetg])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1st_philmetl_w')] += np.sum(weight[one_l_0isrb_ge1st_philmetl])
                hist[sample][tree_name][reference.index('1_el_0isrb_ge1st_w')] += np.sum(weight[one_el_0isrb_ge1st])
                hist[sample][tree_name][reference.index('1_mu_0isrb_ge1st_w')] += np.sum(weight[one_mu_0isrb_ge1st])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt12_r8_w')] += np.sum(weight[one_l_0isrb_ge1sv_pt12_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt12_r9_w')] += np.sum(weight[one_l_0isrb_ge1sv_pt12_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt12_r95_w')] += np.sum(weight[one_l_0isrb_ge1sv_pt12_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt10_r8_w')] += np.sum(weight[one_l_0isrb_ge1sv_pt10_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt10_r9_w')] += np.sum(weight[one_l_0isrb_ge1sv_pt10_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt10_r95_w')] += np.sum(weight[one_l_0isrb_ge1sv_pt10_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt7_r8_w')] += np.sum(weight[one_l_0isrb_ge1sv_pt7_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt7_r9_w')] += np.sum(weight[one_l_0isrb_ge1sv_pt7_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_pt7_r95_w')] += np.sum(weight[one_l_0isrb_ge1sv_pt7_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_dphlep_r8_w')] += np.sum(weight[one_l_0isrb_ge1sv_dphlep_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_dphlep_r9_w')] += np.sum(weight[one_l_0isrb_ge1sv_dphlep_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_dphlep_r95_w')] += np.sum(weight[one_l_0isrb_ge1sv_dphlep_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_dphlep_pt7_r8_w')] += np.sum(weight[one_l_0isrb_ge1sv_dphlep_pt7_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_dphlep_pt7_r9_w')] += np.sum(weight[one_l_0isrb_ge1sv_dphlep_pt7_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_dphlep_pt7_r95_w')] += np.sum(weight[one_l_0isrb_ge1sv_dphlep_pt7_r95])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_chp_r8_w')] += np.sum(weight[one_l_0isrb_ge1sv_chp_r8])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_chp_r9_w')] += np.sum(weight[one_l_0isrb_ge1sv_chp_r9])
                hist[sample][tree_name][reference.index('1l_0isrb_ge1sv_chp_r95_w')] += np.sum(weight[one_l_0isrb_ge1sv_chp_r95])


                print 'finished filling'
    return hist, reference


             

if __name__ == "__main__":

    parser = arg.ArgumentParser(description='receiving file lists for batch jobs')
    parser.add_argument('-f', '--file_list', type=str)
    parser.add_argument('-o', '--out_file', type=str) 

    args = parser.parse_args()

    sample_file = args.file_list
    file_name = args.out_file

    sample_list = pickle.load( open(sample_file, "rb"))

    variables = ['MET', 'MET_phi', 'Eta_SV', 'index_SV_S', 'index_SV_ISR', 'Phi_lep', 'PT_lep', 'Charge_lep', 'ID_lep', 'PDGID_lep', 'index_lep_ISR', 'index_lep_S', 'RISR', 'PTISR', 'SIP3D_lep','MiniIso_lep','MS', 'weight', 'Njet_ISR', 'Njet_S', "Nbjet_ISR", "Nbjet_S", "NSV_ISR", 'NSV_S', 'Nlep', 'Nlep_S', 'Nlep_ISR','PLa', 'PLb', 'PS', 'dphi_lep_S']

    start_b = time.time()    

    sample_arrays, reference_array = make_tables(sample_list, variables, None)
   
    with open(file_name, "wb") as f: 
        pickle.dump(sample_arrays,  f)
        pickle.dump(reference_array, f)
        f.close()
   
    stop_b = time.time()

    print "total: ", stop_b - start_b
 
    print 'finished writing'
