#!/usr/bin/env python

import ROOT as rt
from DataFormats.FWLite import Events, Handle
from rootpy.vector import LorentzVector
from rootpy.tree import Tree, TreeModel, FloatArrayCol, IntCol, FloatCol
from heapq import nlargest

class Event(TreeModel):
    nevents = IntCol()
    npv = IntCol()
    n_AK8 = IntCol()
    n_Higgs = IntCol()
    n_Top = IntCol()
    n_AK8MatchedToHiggs = IntCol()
    n_AK8MatchedToTop = IntCol()
    m_tH = IntCol()
    m_W = FloatArrayCol(10, length_name = "n_Top")
    m_Top = FloatArrayCol(10, length_name = "n_Top")
    m_Higgs = FloatArrayCol(10, length_name = "n_Higgs")

    pt_AK8 = FloatArrayCol(10, length_name='n_AK8')
    eta_AK8 = FloatArrayCol(10, length_name='n_AK8')
    phi_AK8 = FloatArrayCol(10, length_name='n_AK8')
    e_AK8 = FloatArrayCol(10, length_name='n_AK8')

    pt_Higgs = FloatArrayCol(10, length_name='n_Higgs')
    eta_Higgs = FloatArrayCol(10, length_name='n_Higgs')
    phi_Higgs = FloatArrayCol(10, length_name='n_Higgs')
    e_Higgs = FloatArrayCol(10, length_name='n_Higgs')

    pt_Top = FloatArrayCol(10, length_name='n_Top')
    eta_Top = FloatArrayCol(10, length_name='n_Top')
    phi_Top = FloatArrayCol(10, length_name='n_Top')
    e_Top = FloatArrayCol(10, length_name='n_Top')

    pt_W = FloatArrayCol(10, length_name='n_Top')
    eta_W = FloatArrayCol(10, length_name='n_Top')
    phi_W = FloatArrayCol(10, length_name='n_Top')
    e_W = FloatArrayCol(10, length_name='n_Top')

    pt_MatchedHiggs = FloatArrayCol(10, length_name='n_AK8MatchedToHiggs')
    eta_MatchedHiggs = FloatArrayCol(10, length_name='n_AK8MatchedToHiggs')
    phi_MatchedHiggs = FloatArrayCol(10, length_name='n_AK8MatchedToHiggs')
    e_MatchedHiggs = FloatArrayCol(10, length_name='n_AK8MatchedToHiggs')
    m_MatchedHiggs = FloatArrayCol(10, length_name='n_AK8MatchedToHiggs')

    pt_MatchedTop = FloatArrayCol(10, length_name='n_AK8MatchedToTop')
    eta_MatchedTop = FloatArrayCol(10, length_name='n_AK8MatchedToTop')
    phi_MatchedTop = FloatArrayCol(10, length_name='n_AK8MatchedToTop')
    e_MatchedTop = FloatArrayCol(10, length_name='n_AK8MatchedToTop')
    m_MatchedTop = FloatArrayCol(10, length_name='n_AK8MatchedToTop')

    pt_AK8MatchedToHiggs = FloatArrayCol(10, length_name='n_AK8MatchedToHiggs')
    eta_AK8MatchedToHiggs = FloatArrayCol(10, length_name='n_AK8MatchedToHiggs')
    phi_AK8MatchedToHiggs = FloatArrayCol(10, length_name='n_AK8MatchedToHiggs')
    e_AK8MatchedToHiggs = FloatArrayCol(10, length_name='n_AK8MatchedToHiggs')
    m_AK8MatchedToHiggs = FloatArrayCol(10, length_name='n_AK8MatchedToHiggs')

    pt_AK8MatchedToTop = FloatArrayCol(10, length_name='n_AK8MatchedToTop')
    eta_AK8MatchedToTop = FloatArrayCol(10, length_name='n_AK8MatchedToTop')
    phi_AK8MatchedToTop = FloatArrayCol(10, length_name='n_AK8MatchedToTop')
    e_AK8MatchedToTop = FloatArrayCol(10, length_name='n_AK8MatchedToTop')
    m_AK8MatchedToTop = FloatArrayCol(10, length_name='n_AK8MatchedToTop')

TbtH_1200_LH = [
'~/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/TprimeBToTH_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2/TprimeBToTH_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161021_200929/0000/B2GEDMNtuple_1.root',
'~/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/TprimeBToTH_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2/TprimeBToTH_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161021_200929/0000/B2GEDMNtuple_2.root',
'~/eos/cms/store/group/phys_b2g/B2GAnaFW_80X_V2p1/TprimeBToTH_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2/TprimeBToTH_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1_B2GAnaFW_80X_V2p1/161021_200929/0000/B2GEDMNtuple_3.root',
               ]

def Match(inFiles,sample):

    h_npv = Handle("std::int")
    h_genparticles = Handle("vector<reco::GenParticle>")

    h_ak8_pt = Handle("std::vector<float>")
    h_ak8_eta = Handle("std::vector<float>")
    h_ak8_phi = Handle("std::vector<float>")
    h_ak8_e = Handle("std::vector<float>")

    l_npv = ("vertexInfo", "npv")

    l_genparticles  = ("filteredPrunedGenParticles")

    l_ak8_pt = ("jetsAK8CHS", "jetAK8CHSPt")
    l_ak8_eta = ("jetsAK8CHS", "jetAK8CHSEta")
    l_ak8_phi = ("jetsAK8CHS", "jetAK8CHSPhi")
    l_ak8_e  = ("jetsAK8CHS", "jetAK8CHSE")

    fileOut = rt.TFile.Open("jetMatchingTree_Hbb_t_"+sample+"_.root", "RECREATE")
    
    tree = Tree("jetMatchingTree", "Tree for matching AK8 jets to gen H and t", model = Event)
    evts_processed = 0
    nOverallTop = 0
    nOverallHiggs = 0
    nPair = 0
    for f in inFiles:

        events = Events(f) 

        nevents = 0

        for event in events:
            if evts_processed%1000 == 0:
                print 'Events processed: ', evts_processed
            #if int(evts_processed) == 10000: break
            evts_processed += 1
            event.getByLabel(l_npv, h_npv)
            event.getByLabel(l_genparticles, h_genparticles)

            event.getByLabel(l_ak8_pt, h_ak8_pt)
            event.getByLabel(l_ak8_eta, h_ak8_eta)
            event.getByLabel(l_ak8_phi, h_ak8_phi)
            event.getByLabel(l_ak8_e, h_ak8_e)

            if len(h_npv.product()) < 1: continue

            #if len(h_ak8_pt.product()) < 1: continue

            tree.npv = h_npv.product()[0]

            n_AK8 = 0
            n_Higgs = 0
            n_Top = 0
            
            n_MatchedToHiggs = 0
            n_MatchedToTop = 0
            n_AK8MatchedToHiggs = 0
            n_AK8MatchedToTop = 0
            m_tH = 0
            p4_Higgs = []
            p4_Hdau0 = []
            p4_Hdau1 = []
            p4_Top = []
            p4_W = []
            p4_Tdau0 = []
            p4_Tdau1 = []
            p4_AK8 = []
            if h_genparticles.product().size() > 0:
                for igen in xrange(h_genparticles.product().size()):
                        
                    pdgid = h_genparticles.product()[igen].pdgId()
          

                    nDaughters = h_genparticles.product()[igen].numberOfDaughters()
                                
                    if abs(pdgid) == 25:
                        if nDaughters != 2: continue

                        firstDaughter = h_genparticles.product()[igen].daughter(0).pdgId()
                        secondDaughter = h_genparticles.product()[igen].daughter(1).pdgId()
                        if abs(firstDaughter) != 5 or abs(secondDaughter) != 5: continue

                        genPt = h_genparticles.product()[igen].pt() 
                        genEta = h_genparticles.product()[igen].eta()
                        genPhi = h_genparticles.product()[igen].phi()
                        genE = h_genparticles.product()[igen].energy()

                        dau0Pt = h_genparticles.product()[igen].daughter(0).pt()
                        dau0Eta = h_genparticles.product()[igen].daughter(0).eta()
                        dau0Phi = h_genparticles.product()[igen].daughter(0).phi()
                        dau0E = h_genparticles.product()[igen].daughter(0).energy()

                        dau1Pt = h_genparticles.product()[igen].daughter(1).pt()
                        dau1Eta = h_genparticles.product()[igen].daughter(1).eta()
                        dau1Phi = h_genparticles.product()[igen].daughter(1).phi()
                        dau1E = h_genparticles.product()[igen].daughter(1).energy()

                        p4GenHiggs = rt.TLorentzVector()
                        p4GenHiggs.SetPtEtaPhiE(genPt, genEta, genPhi, genE)

                        p4Hdau0 = rt.TLorentzVector()
                        p4Hdau0.SetPtEtaPhiE(dau0Pt, dau0Eta, dau0Phi, dau0E)

                        p4Hdau1 = rt.TLorentzVector()
                        p4Hdau1.SetPtEtaPhiE(dau1Pt, dau1Eta, dau1Phi, dau1E)

                        p4_Higgs.append(p4GenHiggs)
                        p4_Hdau0.append(p4Hdau0)
                        p4_Hdau1.append(p4Hdau1)

                        tree.pt_Higgs[n_Higgs] = p4GenHiggs.Pt()
                        tree.eta_Higgs[n_Higgs] = p4GenHiggs.Eta()
                        tree.phi_Higgs[n_Higgs] = p4GenHiggs.Phi()
                        tree.e_Higgs[n_Higgs] = p4GenHiggs.Energy()
                        tree.m_Higgs[n_Higgs] = p4GenHiggs.Mag()
                      
                        n_Higgs += 1

                    elif abs(pdgid) == 6:
                        if nDaughters != 2: continue
                        
                        firstDaughter = h_genparticles.product()[igen].daughter(0).pdgId()
                        secondDaughter = h_genparticles.product()[igen].daughter(1).pdgId()
                        hadronic = False
                        if abs(firstDaughter) == 24 and abs(secondDaughter) == 5 or abs(firstDaughter) == 5 and abs(secondDaughter) == 24:
                            for idaughter in xrange(nDaughters):
                                if abs(h_genparticles.product()[igen].daughter(idaughter).pdgId()) != 24: continue
                                nGDaughters = h_genparticles.product()[igen].daughter(idaughter).numberOfDaughters()
                                for igDaughter in xrange(nGDaughters):
                                    gDaughter = h_genparticles.product()[igen].daughter(idaughter).daughter(igDaughter).pdgId()
                                    if abs(gDaughter) == 24:
                                        nggDaughters = h_genparticles.product()[igen].daughter(idaughter).daughter(igDaughter).numberOfDaughters()
                                        for iggDaughter in xrange(nggDaughters):
                                            ggDaughter = h_genparticles.product()[igen].daughter(idaughter).daughter(igDaughter).daughter(iggDaughter).pdgId()
                                            if abs(ggDaughter) in xrange(1,7): hadronic = True  
                                            else: hadronic = False                                        
                                    else: 
                                        if abs(gDaughter) in xrange(1,7): hadronic = True
                                        else: hadronic = False
                            if not hadronic: continue

                            genPt = h_genparticles.product()[igen].pt() 
                            genEta = h_genparticles.product()[igen].eta()
                            genPhi = h_genparticles.product()[igen].phi()
                            genE = h_genparticles.product()[igen].energy()

                            dau0Pt = h_genparticles.product()[igen].daughter(0).pt()
                            dau0Eta = h_genparticles.product()[igen].daughter(0).eta()
                            dau0Phi = h_genparticles.product()[igen].daughter(0).phi()
                            dau0E = h_genparticles.product()[igen].daughter(0).energy()

                            dau1Pt = h_genparticles.product()[igen].daughter(1).pt()
                            dau1Eta = h_genparticles.product()[igen].daughter(1).eta()
                            dau1Phi = h_genparticles.product()[igen].daughter(1).phi()
                            dau1E = h_genparticles.product()[igen].daughter(1).energy()

                            p4GenTop = rt.TLorentzVector()
                            p4GenTop.SetPtEtaPhiE(genPt, genEta, genPhi, genE)

                            p4Tdau0 = rt.TLorentzVector()
                            p4Tdau0.SetPtEtaPhiE(dau0Pt, dau0Eta, dau0Phi, dau0E)

                            p4Tdau1 = rt.TLorentzVector()
                            p4Tdau1.SetPtEtaPhiE(dau1Pt, dau1Eta, dau1Phi, dau1E)

                            p4_Top.append(p4GenTop)
                            p4_Tdau0.append(p4Tdau0)
                            p4_Tdau1.append(p4Tdau1)

                            tree.pt_Top[n_Top] = p4GenTop.Pt()
                            tree.eta_Top[n_Top] = p4GenTop.Eta()
                            tree.phi_Top[n_Top] = p4GenTop.Phi()
                            tree.e_Top[n_Top] = p4GenTop.Energy()
                            tree.m_Top[n_Top] = p4GenTop.Mag()
                            
                            if firstDaughter == 24:
                                tree.pt_W[n_Top] = h_genparticles.product()[igen].daughter(0).pt()
                                tree.eta_W[n_Top] = h_genparticles.product()[igen].daughter(0).eta()
                                tree.phi_W[n_Top] = h_genparticles.product()[igen].daughter(0).phi()
                                tree.e_W[n_Top] = h_genparticles.product()[igen].daughter(0).energy()

                                tree.m_W[n_Top] = p4Tdau0.Mag()
                            elif secondDaughter == 24:
                                p4_W.append(p4Tdau1)
                                tree.pt_W[n_Top] = h_genparticles.product()[igen].daughter(1).pt()
                                tree.eta_W[n_Top] = h_genparticles.product()[igen].daughter(1).eta()
                                tree.phi_W[n_Top] = h_genparticles.product()[igen].daughter(1).phi()
                                tree.e_W[n_Top] = h_genparticles.product()[igen].daughter(1).energy()
                               
                                tree.m_W[n_Top] = p4Tdau1.Mag()
                            
                            n_Top += 1

            for ijet in xrange(len(h_ak8_pt.product())):
                jetPt = h_ak8_pt.product()[ijet]
                jetEta = h_ak8_eta.product()[ijet]
                jetPhi = h_ak8_phi.product()[ijet]
                jetE = h_ak8_e.product()[ijet]
                
                if abs(jetPt) < 300: continue
                if abs(jetEta) > 2.4: continue

                p4AK8 = rt.TLorentzVector()
                p4AK8.SetPtEtaPhiE(jetPt, jetEta, jetPhi, jetE)

                p4_AK8.append(p4AK8)

                tree.pt_AK8[n_AK8] = p4AK8.Pt()
                tree.eta_AK8[n_AK8] = p4AK8.Eta()
                tree.phi_AK8[n_AK8] = p4AK8.Phi()
                tree.e_AK8[n_AK8] = p4AK8.Energy()
                n_AK8 += 1
                if p4_Higgs:
                    for iHiggs in xrange(n_Higgs):
                        if p4AK8.DeltaR(p4_Higgs[iHiggs]) > 0.3: continue
                        if p4AK8.DeltaR(p4_Hdau0[iHiggs]) > 0.8: continue
                        if p4AK8.DeltaR(p4_Hdau1[iHiggs]) > 0.8: continue

                        tree.pt_MatchedHiggs[n_AK8MatchedToHiggs] = p4_Higgs[iHiggs].Pt() 
                        tree.eta_MatchedHiggs[n_AK8MatchedToHiggs] = p4_Higgs[iHiggs].Eta() 
                        tree.phi_MatchedHiggs[n_AK8MatchedToHiggs] = p4_Higgs[iHiggs].Phi() 
                        tree.e_MatchedHiggs[n_AK8MatchedToHiggs] = p4_Higgs[iHiggs].Energy() 
                        tree.m_MatchedHiggs[n_AK8MatchedToHiggs] = p4_Higgs[iHiggs].Mag() 
    
                        tree.pt_AK8MatchedToHiggs[n_AK8MatchedToHiggs] = p4AK8.Pt()
                        tree.eta_AK8MatchedToHiggs[n_AK8MatchedToHiggs] = p4AK8.Eta()
                        tree.phi_AK8MatchedToHiggs[n_AK8MatchedToHiggs] = p4AK8.Phi()
                        tree.e_AK8MatchedToHiggs[n_AK8MatchedToHiggs] = p4AK8.Energy()
                        tree.m_AK8MatchedToHiggs[n_AK8MatchedToHiggs] = p4AK8.Mag()
    
                        n_AK8MatchedToHiggs = n_AK8MatchedToHiggs + 1
                if p4_Top:
                    for iTop in xrange(n_Top):
                        if p4AK8.DeltaR(p4_Top[iTop]) > 0.3: continue
                        if p4AK8.DeltaR(p4_Tdau0[iTop]) > 0.8: continue
                        if p4AK8.DeltaR(p4_Tdau1[iTop]) > 0.8: continue
    
                        tree.pt_MatchedTop[n_MatchedToTop] = p4_Top[iTop].Pt() 
                        tree.eta_MatchedTop[n_MatchedToTop] = p4_Top[iTop].Eta() 
                        tree.phi_MatchedTop[n_MatchedToTop] = p4_Top[iTop].Phi() 
                        tree.e_MatchedTop[n_MatchedToTop] = p4_Top[iTop].Energy() 
                        tree.m_MatchedTop[n_MatchedToTop] = p4_Top[iTop].Mag() 
    
                        tree.pt_AK8MatchedToTop[n_AK8MatchedToTop] = p4AK8.Pt()
                        tree.eta_AK8MatchedToTop[n_AK8MatchedToTop] = p4AK8.Eta()
                        tree.phi_AK8MatchedToTop[n_AK8MatchedToTop] = p4AK8.Phi()
                        tree.e_AK8MatchedToTop[n_AK8MatchedToTop] = p4AK8.Energy()
                        tree.m_AK8MatchedToTop[n_AK8MatchedToTop] = p4AK8.Mag()
     
                        n_AK8MatchedToTop = n_AK8MatchedToTop + 1

            if p4_Higgs and p4_Top:
                m_tH = (p4_Higgs[0] + p4_Top[0]).Mag()
                tree.m_tH = m_tH
            if n_AK8MatchedToTop > 0:
                nOverallTop += 1
            if n_AK8MatchedToHiggs > 0:
                nOverallHiggs += 1
            if n_AK8MatchedToTop > 0 and n_AK8MatchedToHiggs > 0:
                nPair += 1
            tree.n_AK8 = n_AK8
            tree.n_Higgs = n_Higgs
            tree.n_Top = n_Top
            tree.n_AK8MatchedToHiggs = n_AK8MatchedToHiggs
            tree.n_AK8MatchedToTop = n_AK8MatchedToTop
            tree.nevents = nevents
            tree.fill()
 
            nevents += 1
    print nOverallTop
    print nOverallHiggs
    print nPair
    tree.Write()
    fileOut.Close()

if __name__ == "__main__":
   Match(TbtH_1200_LH,"TbtH_1200_LH") 
