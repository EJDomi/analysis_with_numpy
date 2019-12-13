import ROOT as rt
from time import strftime, localtime
from collections import OrderedDict
from plotting_susy_cff import plot_configurables as pc
from plotting_susy_cff import sample_configurables as sc
import imp, os
from file_table_functions import *

date = strftime('%d%b%y', localtime())


rt.gROOT.SetBatch()
rt.gROOT.SetStyle('Plain')
rt.gStyle.SetOptTitle(0)
rt.gStyle.SetOptStat(0000)
rt.gStyle.SetOptFit(0111)
rt.gStyle.SetPalette(rt.kBlueRedYellow)
rt.TH1.AddDirectory(rt.kFALSE)

helper   = imp.load_source('fix'     , './help.py')
tdrstyle = imp.load_source('tdrstyle', './tdrstyle.py')
CMS_lumi = imp.load_source('CMS_lumi', './CMS_lumi.py') 

tdrstyle.setTDRStyle()


def make_me_a_canvas():
   can = rt.TCanvas('canvas', 'canvas', 800, 600)
   can.SetLeftMargin(0.15)
   can.SetRightMargin(0.18)
   can.SetBottomMargin(0.15)
   can.SetGridx()
   can.SetGridy()
   can.SetLogz()
   return can


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


def make_new_hists(hists_, hists_2l_):
    temp_new = OrderedDict()
    for sample in hists_:
        temp_new[sample] = OrderedDict()
        for tree in hists_[sample]:
            if '4bd' in sample and sample.split('_')[1] not in tree: continue
            temp_new[sample][tree] = OrderedDict()
            for hist_name, hist in hists_[sample][tree].items():
                new_hist = None
                if 'NS_' in hist_name:
                    new_hist = hist_name.replace('S_', 'S_frac_')
                    temp_new[sample][tree][new_hist] = hist.Clone(hist.GetName().replace('S_','S_frac_'))
                    temp_new[sample][tree][hist_name] = hist.Clone()
                    temp_new[sample][tree][hist_name.replace('S_','Sscale_')] = hist.Clone(hist.GetName().replace('S_','Sscale_'))

                    zero_value = temp_new[sample][tree][new_hist].GetBinContent(1)
                    temp_new[sample][tree][new_hist].Scale(1./zero_value)
                    temp_new[sample][tree][hist_name.replace('S_','Sscale_')].Scale(1./zero_value)

                    temp_inclusive = hists_[sample][tree][hist_name.replace('S_','_')].Clone('temp_S')
                    temp_new[sample][tree][hist_name.replace('S_','scale_')] = hists_[sample][tree][hist_name.replace('S_','_')].Clone(hist.GetName().replace('S_','scale_'))
                    temp_new[sample][tree][hist_name.replace('S_','_')] = hists_[sample][tree][hist_name.replace('S_','_')].Clone()

                    zero_value_inclusive = temp_inclusive.GetBinContent(1)
                    temp_inclusive.Scale(1./ zero_value_inclusive)
                    temp_new[sample][tree][hist_name.replace('S_','scale_')].Scale(1./zero_value_inclusive)
                    
                    temp_new[sample][tree][new_hist].Divide(temp_inclusive)
#                    for ibin in xrange(temp_new[sample][tree][new_hist].GetNbinsX()):
#                        bin_val = temp_new[sample][tree][new_hist].GetBinContent(ibin)
#                        temp_new[sample][tree][new_hist].SetBinContent(ibin, bin_val / zero_value)
                elif 'NISR_' in hist_name:
                    new_hist = hist_name.replace('ISR_', 'ISR_frac_')
                    temp_new[sample][tree][new_hist] = hist.Clone(hist.GetName().replace('ISR_','ISR_frac_'))
                    temp_new[sample][tree][hist_name] = hist.Clone()
                    temp_new[sample][tree][hist_name.replace('ISR_','ISRscale_')] = hist.Clone(hist.GetName().replace('ISR_','ISRscale_'))

                    zero_value = temp_new[sample][tree][new_hist].GetBinContent(1)
                    temp_new[sample][tree][new_hist].Scale(1./zero_value)
                    temp_new[sample][tree][hist_name.replace('ISR_','ISRscale_')].Scale(1./zero_value)

                    temp_inclusive = hists_[sample][tree][hist_name.replace('ISR_','_')].Clone('temp_ISR')

                    zero_value_inclusive = temp_inclusive.GetBinContent(1)
                    temp_inclusive.Scale(1./ zero_value_inclusive)

                    temp_new[sample][tree][new_hist].Divide(temp_inclusive)
#                    for ibin in xrange(temp_new[sample][tree][new_hist].GetNbinsX()):
#                        bin_val = temp_new[sample][tree][new_hist].GetBinContent(ibin)
#                        temp_new[sample][tree][new_hist].SetBinContent(ibin, bin_val / zero_value)


    for sample in hists_:
        for tree in hists_[sample]:
            if '4bd' in sample and sample.split('_')[1] not in tree: continue
            temp_new[sample][tree]['NISR_objects'] = rt.TH1D('NISR_objects_'+sample+'_'+tree, 'N S Objects', 15,  -0.5, 14.5)
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(1),'0')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(2),'1')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(3),'#geq 2')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(4),' ')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(5),'0')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(6),'1')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(7),'#geq 2')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(8),' ')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(9),'0')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(10),'1')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(11),'#geq 2')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(12),' ')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(13),'0')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(14),'1')
            temp_new[sample][tree]['NISR_objects'].GetXaxis().SetBinLabel(temp_new[sample][tree]['NISR_objects'].GetBin(15),'#geq 2')

            for hist_name, hist in hists_[sample][tree].items():
                
                if 'ISR_SV' not in hist_name and 'ISR_bjets' not in hist_name: continue
                print hist_name
                if 'ISR_SV' in hist_name:
                    bin_val_0 = hist.GetBinContent(1)
                    bin_val_1 = hist.GetBinContent(2)
                    bin_val_2plus = 0.
                    for ibin in xrange(3, hist.GetNbinsX()):
                        bin_val_2plus += hist.GetBinContent(ibin)
                    print hist_name, bin_val_0, bin_val_1, bin_val_2plus
                    temp_new[sample][tree]['NISR_objects'].SetBinContent(temp_new[sample][tree]['NISR_objects'].GetBin(1), bin_val_0 / bin_val_0)
                    
                    temp_new[sample][tree]['NISR_objects'].SetBinContent(temp_new[sample][tree]['NISR_objects'].GetBin(2), bin_val_1 / bin_val_0)
                
                    temp_new[sample][tree]['NISR_objects'].SetBinContent(temp_new[sample][tree]['NISR_objects'].GetBin(3), bin_val_2plus / bin_val_0)
                    
                if 'ISR_bjets' in hist_name:
                    bin_val_0 = hist.GetBinContent(1)
                    bin_val_1 = hist.GetBinContent(2)
                    bin_val_2plus = 0.
                    for ibin in xrange(3, hist.GetNbinsX()):
                        bin_val_2plus += hist.GetBinContent(ibin)
                    print hist_name, bin_val_0, bin_val_1, bin_val_2plus
                    temp_new[sample][tree]['NISR_objects'].SetBinContent(temp_new[sample][tree]['NISR_objects'].GetBin(5), bin_val_0 / bin_val_0)
                    
                    temp_new[sample][tree]['NISR_objects'].SetBinContent(temp_new[sample][tree]['NISR_objects'].GetBin(6), bin_val_1 / bin_val_0)
                
                    temp_new[sample][tree]['NISR_objects'].SetBinContent(temp_new[sample][tree]['NISR_objects'].GetBin(7), bin_val_2plus / bin_val_0)
            
    for sample in hists_2l_:
        for tree in hists_2l_[sample]:
            if '4bd' in sample and sample.split('_')[1] not in tree: continue

            for hist_name, hist in hists_2l_[sample][tree].items():
                if 'ISR_SV' not in hist_name and 'ISR_bjets' not in hist_name: continue
                if 'ISR_SV' in hist_name:
                    bin_val_0 = hist.GetBinContent(1)
                    bin_val_1 = hist.GetBinContent(2)
                    bin_val_2plus = 0.
                    for ibin in xrange(3, hist.GetNbinsX()):
                        bin_val_2plus += hist.GetBinContent(ibin)
                    print hist_name, bin_val_0, bin_val_1, bin_val_2plus
                    temp_new[sample][tree]['NISR_objects'].SetBinContent(temp_new[sample][tree]['NISR_objects'].GetBin(9), bin_val_0 / bin_val_0)
                    
                    temp_new[sample][tree]['NISR_objects'].SetBinContent(temp_new[sample][tree]['NISR_objects'].GetBin(10), bin_val_1 / bin_val_0)
                
                    temp_new[sample][tree]['NISR_objects'].SetBinContent(temp_new[sample][tree]['NISR_objects'].GetBin(11), bin_val_2plus / bin_val_0)
                    
                if 'ISR_bjets' in hist_name:
                    bin_val_0 = hist.GetBinContent(1)
                    bin_val_1 = hist.GetBinContent(2)
                    bin_val_2plus = 0.
                    for ibin in xrange(3, hist.GetNbinsX()):
                        bin_val_2plus += hist.GetBinContent(ibin)
                    print hist_name, bin_val_0, bin_val_1, bin_val_2plus
                    temp_new[sample][tree]['NISR_objects'].SetBinContent(temp_new[sample][tree]['NISR_objects'].GetBin(13), bin_val_0 / bin_val_0)
                    
                    temp_new[sample][tree]['NISR_objects'].SetBinContent(temp_new[sample][tree]['NISR_objects'].GetBin(14), bin_val_1 / bin_val_0)
                
                    temp_new[sample][tree]['NISR_objects'].SetBinContent(temp_new[sample][tree]['NISR_objects'].GetBin(15), bin_val_2plus / bin_val_0)


    return temp_new


def make_stacked_plots(hists_, sig_hists_ = None, print_plots = True, suffix_=''):
    '''
    Makes stacked plots following the samples that are given in the histogram dictionary
    '''
    tdrstyle.setTDRStyle()
    if print_plots:
        if not (os.path.isdir('./plots_'+date)): os.mkdir('./plots_'+date)
    
    n_entries = OrderedDict()
    hists_tmp = OrderedDict()
    if sig_hists_:
        sig_hists_tmp = OrderedDict()
    out_dir = os.path.join('./plots_'+date)
    for sample in hists_:
        for tree in hists_[sample]:
            for hist_name, hist in hists_[sample][tree].items():
                if hist_name not in pc: continue
                if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                if hist.InheritsFrom(rt.TH2.Class()): continue                    
                if 'object' not in hist_name: continue 
                n_entries[hist_name] = OrderedDict()
                hists_tmp[hist_name] = OrderedDict()
    for sample in hists_:
        for tree in hists_[sample]:
            for hist_name, hist in hists_[sample][tree].items():
                if hist_name not in pc: continue
                if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                if hist.InheritsFrom(rt.TH2.Class()): continue                   
                if 'object' not in hist_name: continue 
                n_entries[hist_name][sample] = hist.Integral(0,10000)
                hists_tmp[hist_name][sample] = hist
    if sig_hists_:
        for sample in sig_hists_:
            for tree in sig_hists_[sample]:
                for hist_name, hist in sig_hists_[sample][tree].items():
                    if hist_name not in pc: continue
                    if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                    if hist.InheritsFrom(rt.TH2.Class()): continue                    
                    if 'object' not in hist_name: continue 
                    sig_hists_tmp[hist_name] = OrderedDict()
        for sample in sig_hists_:
            for tree in sig_hists_[sample]:
                for hist_name, hist in sig_hists_[sample][tree].items():
                    if hist_name not in pc: continue
                    if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                    if hist.InheritsFrom(rt.TH2.Class()): continue                    
                    if 'object' not in hist_name: continue 
                    sig_hists_tmp[hist_name][sample] = OrderedDict()
        for sample in sig_hists_:
            for tree in sig_hists_[sample]:
                for hist_name, hist in sig_hists_[sample][tree].items():
                    if hist_name not in pc: continue
                    if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                    if hist.InheritsFrom(rt.TH2.Class()): continue  
                    if not hist.GetEntries() > 0.: continue                 
                    if 'T2bW' in sample and '460' not in tree: continue 
                    if 'object' not in hist_name: continue 
                    sig_hists_tmp[hist_name][sample][tree] = hist
    for hist in n_entries:
        n_entries[hist] = OrderedDict(sorted(n_entries[hist].items(), key=lambda x: x[1]))
        can = make_me_a_canvas()
        can.cd() 
        leg = rt.TLegend(0.45,0.86,0.80,0.92,'','brNDC') 
        leg.SetBorderSize(0)
        leg.SetNColumns(2)
        leg.SetTextSize(0.028)
        leg.SetMargin(0.2)
        #leg.SetFillStyle(0)
        #leg2 = rt.TLegend(0.45,0.52,0.80,0.86,'','brNDC') 
        #leg2.SetBorderSize(0)
        #leg2.SetNColumns(1)
        #leg2.SetTextSize(0.028)
        #leg2.SetMargin(0.2)

        #leg2.SetFillStyle(0)
        for sample in n_entries[hist]:
            if 'N' not in hist: hists_tmp[hist][sample].Rebin(1)
            #hists_tmp[hist][sample].Scale(1./hists_tmp[hist][sample].GetBinContent(1))
            hists_tmp[hist][sample].SetLineColor(sc[sample]['color'])
            hists_tmp[hist][sample].SetLineStyle(sc[sample]['style'])
            if sc[sample]['fill']: hists_tmp[hist][sample].SetFillColor(sc[sample]['fill'])
            if sc[sample]['fill_style']: hists_tmp[hist][sample].SetFillStyle(sc[sample]['fill_style'])
            print hist, sample
            leg.AddEntry(hists_tmp[hist][sample], sc[sample]['legend'], 'fl')
        if sig_hists_:
            for sample in sig_hists_tmp[hist]:
                for itr, tree in enumerate(sig_hists_tmp[hist][sample]):
                    #sig_hists_tmp[hist][sample][tree].Scale(1/sig_hists_tmp[hist][sample][tree].GetBinContent(1))
                    sig_hists_tmp[hist][sample][tree].SetLineColor(sc[sample]['color']+itr)
                    sig_hists_tmp[hist][sample][tree].SetLineStyle(sc[sample]['style'])
                    sig_hists_tmp[hist][sample][tree].SetLineWidth(sc[sample]['width'])
                    if sc[sample]['fill']: sig_hists_tmp[hist][sample][tree].SetFillColor(sc[sample]['fill'])
                    if sc[sample]['fill_style']: sig_hists_tmp[hist][sample][tree].SetFillStyle(sc[sample]['fill_style'])
                    print hist, sample, tree
                    if 'Events' in tree:
                        stop_m = '500'
                        neut_m = sample.split('_')[-1]
                        if 'T2bW' in sample:
                            sample_n = '-'.join(sample.split('_')[:2])
                        else:
                            sample_n = '-'.join(sample.split('_')[:3])
                            
                        #leg2.AddEntry(sig_hists_tmp[hist][sample][tree], '#splitline{'+sc[sample]['legend'] + ';}{#kern[-0.1]{M}#lower[-0.05]{#kern[-0.1]{(}#lower[0.3]{#tilde{#lower[-0.1]{t}}}#kern[0.1]{)}}='+stop_m+', M(#tilde{#chi}#lower[-0.1]{_{1}}#kern[-0.7]{#lower[0.3]{^{0}}})='+neut_m+'}', 'fl')
                    else:
                        stop_m = tree.split('_')[1]
                        neut_m = tree.split('_')[2]
                        sample_n = sample.split('_')[0]                
                        #leg2.AddEntry(sig_hists_tmp[hist][sample][tree], '#splitline{'+sc[sample]['legend'] + ';}{#kern[-0.1]{M}#lower[-0.05]{#kern[-0.1]{(}#lower[0.3]{#tilde{#lower[-0.1]{t}}}#kern[0.1]{)}}='+stop_m+', M(#tilde{#chi}#lower[-0.1]{_{1}}#kern[-0.7]{#lower[0.3]{^{0}}})='+neut_m+'}', 'fl')
        can.cd()
        for isam, sample in enumerate(reversed(hists_tmp[hist])):
            hists_tmp[hist][sample].Draw('histsame')
            hists_tmp[hist][sample].GetXaxis().SetLabelFont(42)
            hists_tmp[hist][sample].GetXaxis().SetLabelSize(0.05)
            hists_tmp[hist][sample].GetYaxis().SetTitle(pc[hist]['ylabel'])
            hists_tmp[hist][sample].GetYaxis().CenterTitle()
            hists_tmp[hist][sample].GetYaxis().SetTitleFont(42)
            hists_tmp[hist][sample].GetYaxis().SetTitleSize(0.048)
            hists_tmp[hist][sample].GetYaxis().SetTitleOffset(1.2)
            hists_tmp[hist][sample].GetYaxis().SetLabelFont(42)
            hists_tmp[hist][sample].GetYaxis().SetLabelSize(0.05)
            can.Update()
            if pc[hist]['xmax'] is not None and isam == 0: 
                xmin = pc[hist]['xmin']
                xmax = pc[hist]['xmax']
                hists_tmp[hist][sample].GetXaxis().SetRangeUser(xmin, xmax) 
                can.Update()
        # if pc[hist]['ymax'] is not None: 
        #     ymin = pc[hist]['ymin']
        #     ymax = pc[hist]['ymax']
        #     stack[hist].GetYaxis().SetRangeUser(ymin, ymax) 
        # stack[hist].SetMinimum(0.00001)
        if sig_hists_:
            for sample in sig_hists_tmp[hist]:
                for tree in sig_hists_tmp[hist][sample]:
                    if 'N' not in hist: sig_hists_tmp[hist][sample][tree].Rebin(1)
                    sig_hists_tmp[hist][sample][tree].Draw('histsame')
        can.Update()
        CMS_lumi.writeExtraText = 1
        CMS_lumi.extraText = 'Simulation'
        CMS_lumi.CMS_lumi(can, 0, 10)
        leg.Draw()
        #leg2.Draw()
        tl = rt.TLatex()
        tl.SetTextSize(0.04)
        tl.SetTextFont(42)
        tl.DrawLatex(2.26, 0.00025, "1 lepton")
        tl.DrawLatex(0.5, 0.0004, "N ISR SV")
        tl.DrawLatex(4.2, 0.0004, "N ISR bjets")
        tl.DrawLatex(8.3, 0.0004, "N ISR SV")
        tl.DrawLatex(12.2, 0.0004, "N ISR bjets")
        tl.DrawLatex(10.14, 0.00025, "2 leptons")
        rt.gPad.RedrawAxis()
        can.Update()
        if print_plots:
            if hists_tmp[hist]:
                can.SetLogy()
                first_hist = next(reversed(hists_tmp[hist]))
                ymax = hists_tmp[hist][first_hist].GetMaximum()
                hists_tmp[hist][first_hist].SetMinimum(0.001)
                hists_tmp[hist][first_hist].SetMinimum(0.001)
                hists_tmp[hist][first_hist].SetMaximum(10*ymax)
                #hists_tmp[hist][first_hist].SetMaximum(500)
                can.Update()
                can.SaveAs(out_dir+'/hstack_log_'+hist+'_'+suffix_+'.root')
                can.SaveAs(out_dir+'/hstack_log_'+hist+'_'+suffix_+'.pdf')
                can.SetLogy(0)
                hists_tmp[hist][first_hist].SetMinimum(0.001)
                #hists_tmp[hist][first_hist].SetMaximum(1.5*ymax)
                hists_tmp[hist][first_hist].SetMaximum(1.3)
                can.Update()
                can.SaveAs(out_dir+'/hstack_'+hist+'_'+suffix_+'.root')
                can.SaveAs(out_dir+'/hstack_'+hist+'_'+suffix_+'.pdf')


def make_plots(hists_, sig_hists_ = None, print_plots = True, suffix_=''):
    '''
    Makes plots wihtout stacking the hists, following the samples that are given in the histogram dictionary
    '''
    tdrstyle.setTDRStyle()
    if print_plots:
        if not (os.path.isdir('./plots_'+date)): os.mkdir('./plots_'+date)
    
    n_entries = OrderedDict()
    hists_tmp = OrderedDict()
    if sig_hists_:
        sig_hists_tmp = OrderedDict()
    out_dir = os.path.join('./plots_'+date)
    for sample in hists_:
        for tree in hists_[sample]:
            for hist_name, hist in hists_[sample][tree].items():
                if hist_name not in pc: continue
                if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                if hist.InheritsFrom(rt.TH2.Class()): continue                    
                if 'frac' not in hist_name: continue 
                n_entries[hist_name] = OrderedDict()
                hists_tmp[hist_name] = OrderedDict()
    for sample in hists_:
        for tree in hists_[sample]:
            for hist_name, hist in hists_[sample][tree].items():
                if hist_name not in pc: continue
                if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                if hist.InheritsFrom(rt.TH2.Class()): continue          
                if 'frac' not in hist_name: continue 
                nevts = hist.Integral(0,1000)
                if nevts == 0: continue        
                n_entries[hist_name][sample] = hist.Integral(0,10000)
                hists_tmp[hist_name][sample] = hist
    if sig_hists_:
        for sample in sig_hists_:
            for tree in sig_hists_[sample]:
                for hist_name, hist in sig_hists_[sample][tree].items():
                    if hist_name not in pc: continue
                    if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                    if hist.InheritsFrom(rt.TH2.Class()): continue                    
                    if 'frac' not in hist_name: continue
                    sig_hists_tmp[hist_name] = OrderedDict()
        for sample in sig_hists_:
            for tree in sig_hists_[sample]:
                for hist_name, hist in sig_hists_[sample][tree].items():
                    if hist_name not in pc: continue
                    if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                    if hist.InheritsFrom(rt.TH2.Class()): continue                    
                    if 'frac' not in hist_name: continue
                    sig_hists_tmp[hist_name][sample] = OrderedDict()
        for sample in sig_hists_:
            for tree in sig_hists_[sample]:
                for hist_name, hist in sig_hists_[sample][tree].items():
                    if hist_name not in pc: continue
                    if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                    if hist.InheritsFrom(rt.TH2.Class()): continue  
                    if not hist.GetEntries() > 0.: continue                 
                    if 'T2bW' in sample and '460' not in tree: continue 
                    if 'frac' not in hist_name: continue
                    sig_hists_tmp[hist_name][sample][tree] = hist
    for hist in n_entries:
        n_entries[hist] = OrderedDict(sorted(n_entries[hist].items(), key=lambda x: x[1]))
        can = make_me_a_canvas()
        can.cd() 
        leg = rt.TLegend(0.45,0.86,0.80,0.92,'','brNDC') 
        leg.SetBorderSize(0)
        leg.SetNColumns(2)
        leg.SetTextSize(0.028)
        leg.SetMargin(0.2)
        #leg.SetFillStyle(0)

        leg2 = rt.TLegend(0.45,0.52,0.80,0.86,'','brNDC') 
        leg2.SetBorderSize(0)
        leg2.SetEntrySeparation(0.5)
        leg2.SetNColumns(1)
        leg2.SetTextSize(0.03)
        leg2.SetMargin(0.2)

        leg3 = rt.TLegend(0.45,0.45,0.80,0.52,'','brNDC') 
        leg3.SetBorderSize(0)
        leg3.SetNColumns(1)
        leg3.SetTextSize(0.028)
        leg3.SetMargin(0.2)

        #leg2.SetFillStyle(0)
        for sample in n_entries[hist]:
            if 'N' not in hist: hists_tmp[hist][sample].Rebin(1)
            hists_tmp[hist][sample].SetLineColor(sc[sample]['color'])
            hists_tmp[hist][sample].SetLineStyle(sc[sample]['style'])
            hists_tmp[hist][sample].SetLineWidth(sc[sample]['width'])
            #if sc[sample]['fill']: hists_tmp[hist][sample].SetFillColor(sc[sample]['fill'])
            #if sc[sample]['fill_style']: hists_tmp[hist][sample].SetFillStyle(sc[sample]['fill_style'])
            print hist, sample
            leg.AddEntry(hists_tmp[hist][sample], sc[sample]['legend'], 'fl')
        if sig_hists_:
            for sample in sig_hists_tmp[hist]:
                for itr, tree in enumerate(sig_hists_tmp[hist][sample]):
                    sig_hists_tmp[hist][sample][tree].SetLineColor(sc[sample]['color']+itr)
                    sig_hists_tmp[hist][sample][tree].SetLineStyle(sc[sample]['style'])
                    sig_hists_tmp[hist][sample][tree].SetLineWidth(sc[sample]['width'])
                    if sc[sample]['fill']: sig_hists_tmp[hist][sample][tree].SetFillColor(sc[sample]['fill'])
                    if sc[sample]['fill_style']: sig_hists_tmp[hist][sample][tree].SetFillStyle(sc[sample]['fill_style'])
                    print hist, sample, tree
                    if 'Events' in tree:
                        stop_m = '500'
                        neut_m = sample.split('_')[-1]
                        if 'T2bW' in sample:
                            sample_n = '-'.join(sample.split('_')[:2])
                        else:
                            sample_n = '-'.join(sample.split('_')[:3])
                            
                        leg2.AddEntry(sig_hists_tmp[hist][sample][tree], '#splitline{'+sc[sample]['legend'] + ';}{#kern[-0.1]{M}#lower[-0.05]{#kern[-0.1]{(}#lower[0.3]{#tilde{#lower[-0.1]{t}}}#kern[0.1]{)}}='+stop_m+', M(#tilde{#chi}#lower[-0.1]{_{1}}#kern[-0.7]{#lower[0.3]{^{0}}})='+neut_m+'}', 'fl')
                    else:
                        stop_m = tree.split('_')[1]
                        neut_m = tree.split('_')[2]
                        sample_n = sample.split('_')[0]                
                        leg2.AddEntry(sig_hists_tmp[hist][sample][tree], '#splitline{'+sc[sample]['legend'] + ';}{#kern[-0.1]{M}#lower[-0.05]{#kern[-0.1]{(}#lower[0.3]{#tilde{#lower[-0.1]{t}}}#kern[0.1]{)}}='+stop_m+', M(#tilde{#chi}#lower[-0.1]{_{1}}#kern[-0.7]{#lower[0.3]{^{0}}})='+neut_m+'}', 'fl')
        can.cd()
        for isam, sample in enumerate(reversed(hists_tmp[hist])):
            hists_tmp[hist][sample].Draw('histsame')
            hists_tmp[hist][sample].GetXaxis().SetTitle(pc[hist]['xlabel'])
            hists_tmp[hist][sample].GetXaxis().CenterTitle()
            hists_tmp[hist][sample].GetXaxis().SetTitleFont(42)
            hists_tmp[hist][sample].GetXaxis().SetTitleSize(0.048)
            hists_tmp[hist][sample].GetXaxis().SetTitleOffset(1.15)
            hists_tmp[hist][sample].GetXaxis().SetLabelFont(42)
            hists_tmp[hist][sample].GetXaxis().SetLabelSize(0.05)
            hists_tmp[hist][sample].GetYaxis().SetTitle(pc[hist]['ylabel'])
            hists_tmp[hist][sample].GetYaxis().CenterTitle()
            hists_tmp[hist][sample].GetYaxis().SetTitleFont(42)
            hists_tmp[hist][sample].GetYaxis().SetTitleSize(0.048)
            hists_tmp[hist][sample].GetYaxis().SetTitleOffset(1.2)
            hists_tmp[hist][sample].GetYaxis().SetLabelFont(42)
            hists_tmp[hist][sample].GetYaxis().SetLabelSize(0.05)
            hists_tmp[hist][sample].Draw('histsame')

 
            hists_tmp[hist.replace('S_','ISR_')][sample].SetLineColor(sc[sample]['color'])
            hists_tmp[hist.replace('S_','ISR_')][sample].SetLineStyle(10)
            hists_tmp[hist.replace('S_','ISR_')][sample].SetLineWidth(3)
           
            hists_tmp[hist.replace('S_','ISR_')][sample].Draw('histsame')
 
            
            can.Update()
            if pc[hist]['xmax'] is not None and isam == 0: 
                xmin = pc[hist]['xmin']
                xmax = pc[hist]['xmax']
                hists_tmp[hist][sample].GetXaxis().SetRangeUser(xmin, xmax) 
                can.Update()
        # if pc[hist]['ymax'] is not None: 
        #     ymin = pc[hist]['ymin']
        #     ymax = pc[hist]['ymax']
        #     stack[hist].GetYaxis().SetRangeUser(ymin, ymax) 
        # stack[hist].SetMinimum(0.00001)
        if sig_hists_:
            for sample in sig_hists_tmp[hist]:
                for tree in sig_hists_tmp[hist][sample]:
                    if 'N' not in hist: sig_hists_tmp[hist][sample][tree].Rebin(1)
                    sig_hists_tmp[hist][sample][tree].Draw('histsame')
                    sig_hists_tmp[hist.replace('S_','ISR_')][sample][tree].SetLineColor(sc[sample]['color'])
                    sig_hists_tmp[hist.replace('S_','ISR_')][sample][tree].SetLineStyle(10)
                    sig_hists_tmp[hist.replace('S_','ISR_')][sample][tree].SetLineWidth(2)
           
                    sig_hists_tmp[hist.replace('S_','ISR_')][sample][tree].Draw('histsame')
      
        leg_line_S = rt.TLine(0,0,0,0)
        leg_line_S.SetLineStyle(1)
        leg_line_S.SetLineWidth(3)
        leg_line_S.SetLineColor(rt.kBlack)

        leg_line_ISR = rt.TLine(0,0,0,0)
        leg_line_ISR.SetLineStyle(10)
        leg_line_ISR.SetLineWidth(3)
        leg_line_ISR.SetLineColor(rt.kBlack)
        
        leg3.AddEntry(leg_line_S, 'S System', 'l')
        leg3.AddEntry(leg_line_ISR, 'ISR System', 'l')

        can.Update()
        CMS_lumi.writeExtraText = 1
        CMS_lumi.extraText = 'Simulation'
        CMS_lumi.CMS_lumi(can, 0, 10)
        leg.Draw()
        leg2.Draw()
        leg3.Draw()
        can.Update()
        rt.gPad.RedrawAxis()
        if print_plots:
            if hists_tmp[hist]:
                can.SetLogy()
                first_hist = next(reversed(hists_tmp[hist]))
                ymax = hists_tmp[hist][first_hist].GetMaximum()
                hists_tmp[hist][first_hist].SetMinimum(0.001)
                hists_tmp[hist][first_hist].SetMinimum(0.001)
                #hists_tmp[hist][first_hist].SetMaximum(250*ymax)
                hists_tmp[hist][first_hist].SetMaximum(500)
                can.Update()
                can.SaveAs(out_dir+'/hstack_log_'+hist+'_'+suffix_+'.root')
                can.SaveAs(out_dir+'/hstack_log_'+hist+'_'+suffix_+'.pdf')
                can.SetLogy(0)
                hists_tmp[hist][first_hist].SetMinimum(0.001)
                #hists_tmp[hist][first_hist].SetMaximum(1.5*ymax)
                hists_tmp[hist][first_hist].SetMaximum(1.2)
                can.Update()
                can.SaveAs(out_dir+'/hstack_'+hist+'_'+suffix_+'.root')
                can.SaveAs(out_dir+'/hstack_'+hist+'_'+suffix_+'.pdf')
            else:
                can.SetLogy()
                first_sig = next(iter(sig_hists_tmp[hist]))
                first_hist = next(iter(sig_hists_tmp[hist][first_sig]))
                ymax = sig_hists_tmp[hist][first_sig][first_hist].GetMaximum()
                sig_hists_tmp[hist][first_sig][first_hist].SetMinimum(0.001)
                sig_hists_tmp[hist][first_sig][first_hist].SetMinimum(0.001)
                #sig_hists_tmp[hist][first_sig][first_hist].SetMaximum(250*ymax)
                sig_hists_tmp[hist][first_sig][first_hist].SetMaximum(1)
                can.Update()
                can.SaveAs(out_dir+'/hstack_log_'+hist+'_'+suffix_+'.root')
                can.SaveAs(out_dir+'/hstack_log_'+hist+'_'+suffix_+'.pdf')
                can.SetLogy(0)
                sig_hists_tmp[hist][first_sig][first_hist].SetMinimum(0.001)
                # sig_hists_tmp[hist][first_sig][first_hist].SetMaximum(1.5*ymax)
                sig_hists_tmp[hist][first_sig][first_hist].SetMaximum(1.2)
                can.Update()
                can.SaveAs(out_dir+'/hstack_'+hist+'_'+suffix_+'.root')
                can.SaveAs(out_dir+'/hstack_'+hist+'_'+suffix_+'.pdf')


if __name__ == "__main__":

    signal_file = './output_signal_risr_0p8_mixed_nsv_02Dec19.root'
    sig_hists = read_in_hists(signal_file)

    signal_file = './output_signal_risr_0p8_mixed_nsv_2l_02Dec19.root'
    sig_hists_2l = read_in_hists(signal_file)

    sig_hists_new = make_new_hists(sig_hists, sig_hists_2l)

    #write_hists_to_file(sig_hists_new, './proposal_sig.root') 

    background_file = './output_background_risr_0p8_mixed_nsv_02Dec19.root'
    b_hists = read_in_hists(background_file)

    background_file = './output_background_risr_0p8_mixed_nsv_2l_02Dec19.root'
    b_hists_2l = read_in_hists(background_file)

    b_hists_new = make_new_hists(b_hists, b_hists_2l)

    #write_hists_to_file(b_hists_new, './proposal_b.root')
 
    suffix = 'sv_4bd'
    make_stacked_plots(b_hists_new, sig_hists_new, True, suffix)
    make_plots(b_hists_new, sig_hists_new, True, suffix)
