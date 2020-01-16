import ROOT as rt
from time import strftime, localtime
from collections import OrderedDict
from plotting_susy_cff import plot_configurables as pc
from plotting_susy_cff import sample_configurables as sc
import imp, os

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
   #can.SetGridx()
   #can.SetGridy()
   can.SetLogz()
   return can

def make_2D_plots(hists_, suffix_):
    tdrstyle.setTDRStyle()
    rt.gROOT.SetBatch()
    rt.gROOT.SetStyle('Plain')
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0000)
    rt.gStyle.SetOptFit(0111)
    rt.gStyle.SetPalette(rt.kBlueRedYellow)
    if not (os.path.isdir('./plots_'+date)): os.mkdir('./plots_'+date)
    for sample in hists_:
        if not (os.path.isdir('./plots_'+date+'/'+sample)): os.mkdir('./plots_'+date+'/'+sample)
        for tree in hists_[sample]:
            if not (os.path.isdir('./plots_'+date+'/'+sample+'/'+tree)): os.mkdir(os.path.join('./plots_'+date, sample, tree))
            out_dir = os.path.join('./plots_'+date, sample, tree)
            for hist_name, hist in hists_[sample][tree].items():
                if not hist.InheritsFrom(rt.TH2.Class()): continue
                if hist_name not in pc: continue
                can = make_me_a_canvas()
                can.cd() 
                if 'N' not in hist_name.split('_')[0]: hist.RebinX(8)
                if 'N' not in hist_name.split('_')[1]: hist.RebinY(8)
                hist.Draw("COLZ")
                hist.GetXaxis().CenterTitle()
                hist.GetXaxis().SetTitleFont(42)
                hist.GetXaxis().SetTitleSize(0.06)
                hist.GetXaxis().SetTitleOffset(1.06)
                hist.GetXaxis().SetLabelFont(42)
                hist.GetXaxis().SetLabelSize(0.05)
                hist.GetXaxis().SetTitle(pc[hist_name]['xlabel'])
                hist.GetYaxis().CenterTitle()
                hist.GetYaxis().SetTitleFont(42)
                hist.GetYaxis().SetTitleSize(0.06)
                hist.GetYaxis().SetTitleOffset(1.12)
                hist.GetYaxis().SetLabelFont(42)
                hist.GetYaxis().SetLabelSize(0.05)
                hist.GetYaxis().SetTitle(pc[hist_name]['ylabel'])
                hist.GetZaxis().CenterTitle()
                hist.GetZaxis().SetTitleFont(42)
                hist.GetZaxis().SetTitleSize(0.06)
                hist.GetZaxis().SetTitleOffset(1.1)
                hist.GetZaxis().SetLabelFont(42)
                hist.GetZaxis().SetLabelSize(0.05)
                hist.GetZaxis().SetTitle("a. u.")
                if pc[hist_name]['xmax'] is not None: 
                    xmin = pc[hist_name]['xmin']
                    xmax = pc[hist_name]['xmax']
                    hist.GetXaxis().SetRangeUser(xmin, xmax) 
                if pc[hist_name]['ymax'] is not None: 
                    ymin = pc[hist_name]['ymin']
                    ymax = pc[hist_name]['ymax']
                    hist.GetYaxis().SetRangeUser(ymin, ymax) 
                #hist.GetZaxis().SetRangeUser(0.9*hist.GetMinimum(0.0),1.1*hist.GetMaximum())
                hist.GetZaxis().SetRangeUser(0.001,1.1*hist.GetMaximum())
                CMS_lumi.writeExtraText = 1
                CMS_lumi.extraText = 'Simulation'
                CMS_lumi.CMS_lumi(can, 0, 10)
                l = rt.TLatex()
                l.SetTextFont(42)
                l.SetNDC()
                l.SetTextSize(0.055)
                l.SetTextFont(42)
                if 'SMS' in sample:
                    l.DrawLatex(0.45,0.943,' '.join(sc[sample]['legend'].split('-')[-2:])+' '+tree)
                else:
                    l.DrawLatex(0.6,0.943,sc[sample]['legend'])
                can.SetLogz()
                can.Update()
                can.SaveAs(out_dir+'/h_'+hist.GetName()+'_'+suffix_+'.root')
                can.SaveAs(out_dir+'/h_'+hist.GetName()+'_'+suffix_+'.pdf')


def make_overlay_plot(hists_, suffix_):
    print 'for posterity'
    hists_tmp = OrderedDict()
    if not (os.path.isdir('./plots_'+date)): os.mkdir('./plots_'+date)
    out_dir = os.path.join('./plots_'+date)

    for sample in hists_:
        for tree in hists_[sample]:
            for hist_name, hist in hists_[sample][tree].items():
                if hist_name not in pc: continue
                if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                if hist.InheritsFrom(rt.TH2.Class()): continue                    
                hists_tmp[hist_name] = OrderedDict()

    for sample in hists_:
        for tree in hists_[sample]:
            for hist_name, hist in hists_[sample][tree].items():
                if hist_name not in pc: continue
                if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                if hist.InheritsFrom(rt.TH2.Class()): continue                    
                hists_tmp[hist_name][sample] = OrderedDict()
    for sample in hists_:
        for tree in hists_[sample]:
            for hist_name, hist in hists_[sample][tree].items():
                if hist_name not in pc: continue
                if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                if hist.InheritsFrom(rt.TH2.Class()): continue                    
                hists_tmp[hist_name][sample][tree] = hist
    for hist in hists_tmp:
        can = make_me_a_canvas()
        can.cd() 
        leg = rt.TLegend(0.2,0.73,0.75,0.93,'','brNDC') 
        leg.SetBorderSize(0)
        leg.SetNColumns(2)
        leg.SetTextSize(0.02)
        leg.SetMargin(0.2)
        leg.SetFillStyle(0)
        for sample in hists_tmp[hist]:
            for itr, tree in enumerate(hists_tmp[hist][sample]):
                if int(hists_tmp[hist][sample][tree].GetEntries()) == 0: continue
                hists_tmp[hist][sample][tree].SetLineColor(sc[sample]['color']+itr*2)
                hists_tmp[hist][sample][tree].SetLineStyle(sc[sample]['style'])
                hists_tmp[hist][sample][tree].SetLineWidth(sc[sample]['width'])
                if sc[sample]['fill']: hists_tmp[hist][sample][tree].SetFillColor(sc[sample]['fill'])
                if sc[sample]['fill_style']: hists_tmp[hist][sample][tree].SetFillStyle(sc[sample]['fill_style'])
                print hist, sample, tree
                stop_m = tree.split('_')[1]
                neut_m = tree.split('_')[2]
                if 'SMS' in sample:
                    leg.AddEntry(hists_tmp[hist][sample][tree], sc[sample]['legend']+' M(#tilde{t})='+stop_m+', M(#tilde{#chi}_{1}^{0})='+neut_m, 'fl')
                else:
                    leg.AddEntry(hists_tmp[hist][sample], sc[sample]['legend'], 'fl')

                if 'N' not in hist: hists_tmp[hist][sample][tree].Rebin(1)
                #hists_tmp[hist][sample][tree].Scale(1/hists_tmp[hist][sample][tree].Integral())
   
                hists_tmp[hist][sample][tree].GetXaxis().SetTitle(pc[hist]['xlabel'])
                hists_tmp[hist][sample][tree].GetXaxis().CenterTitle()
                hists_tmp[hist][sample][tree].GetXaxis().SetTitleFont(42)
                hists_tmp[hist][sample][tree].GetXaxis().SetTitleSize(0.06)
                hists_tmp[hist][sample][tree].GetXaxis().SetTitleOffset(1.06)
                hists_tmp[hist][sample][tree].GetXaxis().SetLabelFont(42)
                hists_tmp[hist][sample][tree].GetXaxis().SetLabelSize(0.05)
                hists_tmp[hist][sample][tree].GetYaxis().SetTitle(pc[hist]['ylabel'])
                hists_tmp[hist][sample][tree].GetYaxis().CenterTitle()
                hists_tmp[hist][sample][tree].GetYaxis().SetTitleFont(42)
                hists_tmp[hist][sample][tree].GetYaxis().SetTitleSize(0.06)
                hists_tmp[hist][sample][tree].GetYaxis().SetTitleOffset(1.12)
                hists_tmp[hist][sample][tree].GetYaxis().SetLabelFont(42)
                hists_tmp[hist][sample][tree].GetYaxis().SetLabelSize(0.05)

                ymax = hists_tmp[hist][sample][tree].GetMaximum()
                #hists_tmp[hist][sample][tree].SetMinimum(0.01)
                #hists_tmp[hist][sample][tree].SetMaximum(250*ymax)

                hists_tmp[hist][sample][tree].SetMinimum(0.01)
                hists_tmp[hist][sample][tree].SetMaximum(2.*ymax)

                hists_tmp[hist][sample][tree].Draw('histsame')
        can.cd()
        CMS_lumi.writeExtraText = 1
        CMS_lumi.extraText = 'Simulation'
        CMS_lumi.CMS_lumi(can, 0, 10)
        leg.Draw()
        can.SetLogy()
        can.Update()
        #can.SaveAs(out_dir+'/hoverlay_log_'+hist+'_'+suffix_+'.root')
        #can.SaveAs(out_dir+'/hoverlay_log_'+hist+'_'+suffix_+'.pdf')
        can.SetLogy(0)
        can.Update()
        can.SaveAs(out_dir+'/hoverlay_'+hist+'_'+suffix_+'.root')
        can.SaveAs(out_dir+'/hoverlay_'+hist+'_'+suffix_+'.pdf')

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
                n_entries[hist_name] = OrderedDict()
                hists_tmp[hist_name] = OrderedDict()
    for sample in hists_:
        for tree in hists_[sample]:
            for hist_name, hist in hists_[sample][tree].items():
                if hist_name not in pc: continue
                if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                if hist.InheritsFrom(rt.TH2.Class()): continue           
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
                    sig_hists_tmp[hist_name] = OrderedDict()
        for sample in sig_hists_:
            for tree in sig_hists_[sample]:
                for hist_name, hist in sig_hists_[sample][tree].items():
                    if hist_name not in pc: continue
                    if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                    if hist.InheritsFrom(rt.TH2.Class()): continue                    
                    sig_hists_tmp[hist_name][sample] = OrderedDict()
        for sample in sig_hists_:
            for tree in sig_hists_[sample]:
                for hist_name, hist in sig_hists_[sample][tree].items():
                    if hist_name not in pc: continue
                    if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                    if hist.InheritsFrom(rt.TH2.Class()): continue  
                    if not hist.GetEntries() > 0.: continue                 
                    if 'T2bW' in sample and '460' not in tree: continue 
                    sig_hists_tmp[hist_name][sample][tree] = hist
    for hist in n_entries:
        n_entries[hist] = OrderedDict(sorted(n_entries[hist].items(), key=lambda x: x[1]))
        can = make_me_a_canvas()
        can.cd() 
        leg = rt.TLegend(0.45,0.86,0.8,0.92,'','brNDC') 
        leg.SetBorderSize(0)
        leg.SetNColumns(2)
        leg.SetTextSize(0.024)
        leg.SetMargin(0.2)
        #leg.SetFillStyle(0)
        leg2 = rt.TLegend(0.45,0.58,0.80,0.86,'','brNDC') 
        leg2.SetBorderSize(0)
        leg2.SetEntrySeparation(0.5)
        leg2.SetNColumns(1)
        leg2.SetTextSize(0.024)
        leg2.SetMargin(0.2)
        #leg2.SetFillStyle(0)
        for sample in n_entries[hist]:
            if 'N' not in hist: hists_tmp[hist][sample].Rebin(1)
            hists_tmp[hist][sample].Scale(1./n_entries[hist][sample])
            hists_tmp[hist][sample].SetLineColor(sc[sample]['color'])
            hists_tmp[hist][sample].SetLineStyle(sc[sample]['style'])
            hists_tmp[hist][sample].SetLineWidth(sc[sample]['width'])
            if sc[sample]['fill']: hists_tmp[hist][sample].SetFillColor(sc[sample]['fill'])
            if sc[sample]['fill_style']: hists_tmp[hist][sample].SetFillStyle(sc[sample]['fill_style'])
            print hist, sample
            leg.AddEntry(hists_tmp[hist][sample], sc[sample]['legend'], 'fl')
        if sig_hists_:
            for sample in sig_hists_tmp[hist]:
                for itr, tree in enumerate(sig_hists_tmp[hist][sample]):
                    sig_hists_tmp[hist][sample][tree].Scale(1/sig_hists_tmp[hist][sample][tree].Integral())
                    sig_hists_tmp[hist][sample][tree].SetLineColor(sc[sample]['color']+itr+2)
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
                        leg2.AddEntry(sig_hists_tmp[hist][sample][tree], sc[sample]['legend']+'; #kern[-0.1]{M}#lower[-0.05]{#kern[-0.1]{(}#lower[0.3]{#tilde{#lower[-0.1]{t}}}#kern[0.1]{)}}='+stop_m+', M(#tilde{#chi}#lower[-0.1]{_{1}}#kern[-0.7]{#lower[0.3]{^{0}}})='+neut_m, 'fl')
                    else:
                        stop_m = tree.split('_')[1]
                        neut_m = tree.split('_')[2]
                        sample_n = sample.split('_')[0]                
                        #leg2.AddEntry(sig_hists_tmp[hist][sample][tree], '#splitline{'+sc[sample]['legend'] + ';}{#kern[-0.1]{M}#lower[-0.05]{#kern[-0.1]{(}#lower[0.3]{#tilde{#lower[-0.1]{t}}}#kern[0.1]{)}}='+stop_m+', M(#tilde{#chi}#lower[-0.1]{_{1}}#kern[-0.7]{#lower[0.3]{^{0}}})='+neut_m+'}', 'fl')
                        leg2.AddEntry(sig_hists_tmp[hist][sample][tree], sc[sample]['legend']+'; #kern[-0.1]{M}#lower[-0.05]{#kern[-0.1]{(}#lower[0.3]{#tilde{#lower[-0.1]{t}}}#kern[0.1]{)}}='+stop_m+', M(#tilde{#chi}#lower[-0.1]{_{1}}#kern[-0.7]{#lower[0.3]{^{0}}})='+neut_m, 'fl')
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
        #l = rt.TLine(20., 0., 20., 0.2)
        #l.SetLineColor(13)
        #l.SetLineStyle(9)
        #l.SetLineWidth(1)
        #l.Draw()
        leg.Draw()
        leg2.Draw()
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
                can.SaveAs(out_dir+'/hoverlay_log_'+hist+'_'+suffix_+'.root')
                can.SaveAs(out_dir+'/hoverlay_log_'+hist+'_'+suffix_+'.pdf')
                can.SetLogy(0)
                hists_tmp[hist][first_hist].SetMinimum(0.001)
                hists_tmp[hist][first_hist].SetMaximum(2.0*ymax)
                #hists_tmp[hist][first_hist].SetMaximum(0.3)
                can.Update()
                can.SaveAs(out_dir+'/hoverlay_'+hist+'_'+suffix_+'.root')
                can.SaveAs(out_dir+'/hoverlay_'+hist+'_'+suffix_+'.pdf')
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
                can.SaveAs(out_dir+'/hoverlay_log_'+hist+'_'+suffix_+'.root')
                can.SaveAs(out_dir+'/hoverlay_log_'+hist+'_'+suffix_+'.pdf')
                can.SetLogy(0)
                sig_hists_tmp[hist][first_sig][first_hist].SetMinimum(0.001)
                # sig_hists_tmp[hist][first_sig][first_hist].SetMaximum(1.5*ymax)
                sig_hists_tmp[hist][first_sig][first_hist].SetMaximum(0.3)
                can.Update()
                can.SaveAs(out_dir+'/hoverlay_'+hist+'_'+suffix_+'.root')
                can.SaveAs(out_dir+'/hoverlay_'+hist+'_'+suffix_+'.pdf')

    #return stack, leg

def make_stacked_plots(hists_, sig_hists_ = None, print_plots = True, suffix_=''):
    '''
    Makes stacked plots following the samples that are given in the histogram dictionary
    '''
    tdrstyle.setTDRStyle()
    if print_plots:
        if not (os.path.isdir('./plots_'+date)): os.mkdir('./plots_'+date)
    stack = OrderedDict()
    
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
                stack[hist_name] = rt.THStack('stack','')
                n_entries[hist_name] = OrderedDict()
                hists_tmp[hist_name] = OrderedDict()
    for sample in hists_:
        for tree in hists_[sample]:
            for hist_name, hist in hists_[sample][tree].items():
                if hist_name not in pc: continue
                if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                if hist.InheritsFrom(rt.TH2.Class()): continue                    
                n_entries[hist_name][sample] = hist.Integral(0,10000)
                hists_tmp[hist_name][sample] = hist
    if sig_hists_:
        for sample in sig_hists_:
            for tree in sig_hists_[sample]:
                for hist_name, hist in sig_hists_[sample][tree].items():
                    if hist_name not in pc: continue
                    if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                    if hist.InheritsFrom(rt.TH2.Class()): continue                    
                    sig_hists_tmp[hist_name] = OrderedDict()
        for sample in sig_hists_:
            for tree in sig_hists_[sample]:
                for hist_name, hist in sig_hists_[sample][tree].items():
                    if hist_name not in pc: continue
                    if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                    if hist.InheritsFrom(rt.TH2.Class()): continue                    
                    sig_hists_tmp[hist_name][sample] = OrderedDict()
        for sample in sig_hists_:
            for tree in sig_hists_[sample]:
                for hist_name, hist in sig_hists_[sample][tree].items():
                    if hist_name not in pc: continue
                    if not hist.InheritsFrom(rt.TH1.Class()): continue                    
                    if hist.InheritsFrom(rt.TH2.Class()): continue  
                    if not hist.GetEntries() > 0.: continue                 
                    if 'T2bW' in sample and '460' not in tree: continue 
                    sig_hists_tmp[hist_name][sample][tree] = hist
    for hist in n_entries:
        n_entries[hist] = OrderedDict(sorted(n_entries[hist].items(), key=lambda x: x[1]))
        can = make_me_a_canvas()
        can.cd() 
        leg = rt.TLegend(0.45,0.86,0.8,0.92,'','brNDC') 
        leg.SetBorderSize(0)
        leg.SetNColumns(2)
        leg.SetTextSize(0.024)
        leg.SetMargin(0.2)
        #leg.SetFillStyle(0)
        leg2 = rt.TLegend(0.45,0.72,0.80,0.86,'','brNDC') 
        leg2.SetBorderSize(0)
        leg2.SetNColumns(1)
        leg2.SetTextSize(0.024)
        leg2.SetMargin(0.2)
        #leg2.SetFillStyle(0)
        for sample in n_entries[hist]:
            if 'N' not in hist: hists_tmp[hist][sample].Rebin(1)
            # hists_tmp[hist][sample].Scale(1./n_entries[hist][sample])
            hists_tmp[hist][sample].SetLineColor(sc[sample]['color'])
            hists_tmp[hist][sample].SetLineStyle(sc[sample]['style'])
            if sc[sample]['fill']: hists_tmp[hist][sample].SetFillColor(sc[sample]['fill'])
            if sc[sample]['fill_style']: hists_tmp[hist][sample].SetFillStyle(1001)
            print hist, sample
            stack[hist].Add(hists_tmp[hist][sample])
            leg.AddEntry(hists_tmp[hist][sample], sc[sample]['legend'], 'fl')
        if sig_hists_:
            for sample in sig_hists_tmp[hist]:
                for itr, tree in enumerate(sig_hists_tmp[hist][sample]):
                    # sig_hists_tmp[hist][sample][tree].Scale(1/sig_hists_tmp[hist][sample][tree].Integral())
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
                        leg2.AddEntry(sig_hists_tmp[hist][sample][tree], sc[sample]['legend']+'; #kern[-0.1]{M}#lower[-0.05]{#kern[-0.1]{(}#lower[0.3]{#tilde{#lower[-0.1]{t}}}#kern[0.1]{)}}='+stop_m+', M(#tilde{#chi}#lower[-0.1]{_{1}}#kern[-0.7]{#lower[0.3]{^{0}}})='+neut_m, 'fl')
                    else:
                        stop_m = tree.split('_')[1]
                        neut_m = tree.split('_')[2]
                        sample_n = sample.split('_')[0]                
                        #leg2.AddEntry(sig_hists_tmp[hist][sample][tree], '#splitline{'+sc[sample]['legend'] + ';}{#kern[-0.1]{M}#lower[-0.05]{#kern[-0.1]{(}#lower[0.3]{#tilde{#lower[-0.1]{t}}}#kern[0.1]{)}}='+stop_m+', M(#tilde{#chi}#lower[-0.1]{_{1}}#kern[-0.7]{#lower[0.3]{^{0}}})='+neut_m+'}', 'fl')
                        leg2.AddEntry(sig_hists_tmp[hist][sample][tree], sc[sample]['legend']+'; #kern[-0.1]{M}#lower[-0.05]{#kern[-0.1]{(}#lower[0.3]{#tilde{#lower[-0.1]{t}}}#kern[0.1]{)}}='+stop_m+', M(#tilde{#chi}#lower[-0.1]{_{1}}#kern[-0.7]{#lower[0.3]{^{0}}})='+neut_m, 'fl')
        can.cd()
        stack[hist].Draw('hist')
        stack[hist].GetXaxis().SetTitle(pc[hist]['xlabel'])
        stack[hist].GetXaxis().CenterTitle()
        stack[hist].GetXaxis().SetTitleFont(42)
        stack[hist].GetXaxis().SetTitleSize(0.06)
        stack[hist].GetXaxis().SetTitleOffset(1.06)
        stack[hist].GetXaxis().SetLabelFont(42)
        stack[hist].GetXaxis().SetLabelSize(0.05)
        stack[hist].GetYaxis().SetTitle(pc[hist]['ylabel'])
        stack[hist].GetYaxis().CenterTitle()
        stack[hist].GetYaxis().SetTitleFont(42)
        stack[hist].GetYaxis().SetTitleSize(0.06)
        stack[hist].GetYaxis().SetTitleOffset(1.12)
        stack[hist].GetYaxis().SetLabelFont(42)
        stack[hist].GetYaxis().SetLabelSize(0.05)
        can.Update()
        if pc[hist]['xmax'] is not None: 
            xmin = pc[hist]['xmin']
            xmax = pc[hist]['xmax']
            stack[hist].GetXaxis().SetRangeUser(xmin, xmax) 
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
        leg2.Draw()
        can.Update()
        if print_plots:
            can.SetLogy()
            ymax = stack[hist].GetMaximum()
            stack[hist].SetMinimum(0.001)
            stack[hist].SetMaximum(2000*ymax)
            #stack[hist].SetMaximum(1.0)
            can.Update()
            can.SaveAs(out_dir+'/hstack_log_'+hist+'_'+suffix_+'.root')
            can.SaveAs(out_dir+'/hstack_log_'+hist+'_'+suffix_+'.pdf')
            can.SetLogy(0)
            stack[hist].SetMinimum(0.001)
            stack[hist].SetMaximum(2.0*ymax)
            can.Update()
            can.SaveAs(out_dir+'/hstack_'+hist+'_'+suffix_+'.root')
            can.SaveAs(out_dir+'/hstack_'+hist+'_'+suffix_+'.pdf')

    #return stack, leg

def make_1D_plots(hists_, suffix_):
    print 'some more posterity here'
    tdrstyle.setTDRStyle()
    if not (os.path.isdir('./plots_'+date)): os.mkdir('./plots_'+date)
    for sample in hists_:
        if not (os.path.isdir('./plots_'+date+'/'+sample)): os.mkdir('./plots_'+date+'/'+sample)
        for tree in hists_[sample]:
            if not (os.path.isdir('./plots_'+date+'/'+sample+'/'+tree)): os.mkdir(os.path.join('./plots_'+date, sample, tree))
            out_dir = os.path.join('./plots_'+date, sample, tree)
            for hist_name, hist in hists_[sample][tree].items():
                if not hist.InheritsFrom(rt.TH1.Class()): continue                   
                if hist.InheritsFrom(rt.TH2.Class()): continue 
                if hist_name not in pc: continue
                can = make_me_a_canvas()
                can.cd()
                hist.Draw('hist')
                if 'N' not in hist_name: hist.Rebin(8) 
                hist.SetLineColor(sc[sample]['color'])
                hist.SetLineStyle(sc[sample]['style'])
                hist.SetLineWidth(sc[sample]['width'])
                if sc[sample]['fill']: hist.SetFillColor(sc[sample]['fill'])
                if sc[sample]['fill_style']: hist.SetFillStyle(sc[sample]['fill_style'])
                hist.GetXaxis().CenterTitle()
                hist.GetXaxis().SetTitleFont(42)
                hist.GetXaxis().SetTitleSize(0.06)
                hist.GetXaxis().SetTitleOffset(1.06)
                hist.GetXaxis().SetLabelFont(42)
                hist.GetXaxis().SetLabelSize(0.05)
                hist.GetXaxis().SetTitle(pc[hist_name]['xlabel'])
                hist.GetYaxis().CenterTitle()
                hist.GetYaxis().SetTitleFont(42)
                hist.GetYaxis().SetTitleSize(0.06)
                hist.GetYaxis().SetTitleOffset(1.12)
                hist.GetYaxis().SetLabelFont(42)
                hist.GetYaxis().SetLabelSize(0.05)
                hist.GetYaxis().SetTitle(pc[hist_name]['ylabel'])
                hist.GetYaxis().SetRangeUser(0.001,1.2*hist.GetMaximum())
                CMS_lumi.writeExtraText = 1
                CMS_lumi.extraText = 'Simulation'
                CMS_lumi.CMS_lumi(can, 0, 10)
                l = rt.TLatex()
                l.SetTextFont(42)
                l.SetNDC()
                l.SetTextSize(0.055)
                l.SetTextFont(42)
                l.DrawLatex(0.41,0.943,sc[sample]['legend'])
      
                can.Update()
                can.SaveAs(out_dir+'/h_'+hist.GetName()+'_'+suffix_+'.root')
                can.SaveAs(out_dir+'/h_'+hist.GetName()+'_'+suffix_+'.pdf')

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

def make_new_hists(hists_):
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
                    temp_new[sample][tree][new_hist].Divide(hists_[sample][tree][hist_name.replace('S_','_')])
                    zero_value = temp_new[sample][tree][new_hist].GetBinContent(1)
                    temp_new[sample][tree][new_hist].Scale(1./zero_value)
#                    for ibin in xrange(temp_new[sample][tree][new_hist].GetNbinsX()):
#                        bin_val = temp_new[sample][tree][new_hist].GetBinContent(ibin)
#                        temp_new[sample][tree][new_hist].SetBinContent(ibin, bin_val / zero_value)
                elif 'NISR_' in hist_name:
                    new_hist = hist_name.replace('ISR_', 'ISR_frac_')
                    temp_new[sample][tree][new_hist] = hist.Clone(hist.GetName().replace('ISR_','ISR_frac_'))
                    temp_new[sample][tree][new_hist].Divide(hists_[sample][tree][hist_name.replace('ISR_','_')])
                    zero_value = temp_new[sample][tree][new_hist].GetBinContent(1)
                    temp_new[sample][tree][new_hist].Scale(1./zero_value)
#                    for ibin in xrange(temp_new[sample][tree][new_hist].GetNbinsX()):
#                        bin_val = temp_new[sample][tree][new_hist].GetBinContent(ibin)
#                        temp_new[sample][tree][new_hist].SetBinContent(ibin, bin_val / zero_value)
    hists_.update(temp_new)
    return hists_

def read_in_hists_opt(in_file_):
    in_file = rt.TFile(in_file_, 'r')
    hists = OrderedDict()
    for key in in_file.GetListOfKeys():
        key_name = key.GetName()
        print key_name
        key_split = key_name.split('__')
        sample = key_split[0]
        tree = key_split[1]
        hist_name = key_split[2]
        
        if sample not in hists.keys():
            hists[sample] = OrderedDict()
        if tree not in hists[sample].keys():
            hists[sample][tree] = OrderedDict()

        hists[sample][tree][hist_name] = in_file.Get(key_name)
    return hists


if __name__ == "__main__":
    #signal_file = './output_signal_risr_95_mixed_10Jan20.root'
    signal_file = './output_signal_risr_mixed_10Jan20.root'
    sig_hists = read_in_hists(signal_file)
    #suffix = '0p95'
    suffix = 'all'
    #sig_hists_new = make_new_hists(sig_hists)
    #print sig_hists_new
    make_1D_plots(sig_hists, suffix)

    #background_file = './output_background_risr_95_mixed_10Jan20.root'
    background_file = './output_background_risr_mixed_10Jan20.root'
    b_hists = read_in_hists(background_file)
    suffix = 'all'
    #b_hists = make_new_hists(b_hists)
    make_1D_plots(b_hists, suffix)
    make_stacked_plots(b_hists, sig_hists, True, suffix)
    make_plots(b_hists, sig_hists, True, suffix)
    make_2D_plots(sig_hists, suffix)
    make_2D_plots(b_hists, suffix)
