
import os
import array
import glob
import math
import ROOT
import sys
from ROOT import *
from array import *
from optparse import OptionParser
parser = OptionParser()

parser.add_option('-s', '--set', metavar='F', type='string', action='store',
                  default	=	'data',
                  dest		=	'set',
                  help		=	'data or QCD')

parser.add_option('-c', '--cuts', metavar='F', type='string', action='store',
                  default	=	'rate_default',
                  dest		=	'cuts',
                  help		=	'Cuts type (ie default, rate, etc)')

(options, args) = parser.parse_args()
rootdir="rootfiles/"
import Wprime_Functions	
from Wprime_Functions import *

gROOT.Macro("rootlogon.C")
gROOT.LoadMacro("insertlogo.C+")

BTR = BTR_Init('Bifpoly',options.cuts,'')
BTR_err = BTR_Init('Bifpoly_err',options.cuts,'')
fittitles = ["pol0","pol2","pol3","FIT","expofit"]
fits = []
for fittitle in fittitles:
	fits.append(BTR_Init(fittitle,options.cuts,''))

leg1 = TLegend(0.45,0.57,.84,.78)
leg1.SetFillColor(0)
leg1.SetBorderSize(0)

leg2 = TLegend(0.,0.,1.,1.)
leg2.SetFillColor(0)
leg2.SetBorderSize(0)

#output = ROOT.TFile( "fitting.root", "recreate" )
#output.cd()
c1 = TCanvas('c1', 'Tagrate numerator and deominator', 1000, 1300)
c4 = TCanvas('c4', 'Pt fitted tagrate in 0.0 < Eta <0.5', 800, 500)
c5 = TCanvas('c5', 'Pt fitted tagrate in 0.5 < Eta <1.15', 800, 500)
c6 = TCanvas('c6', 'Pt fitted tagrate in 1.15 < Eta <2.4', 800, 500)
c7 = TCanvas('c7', 'tagged vs signal', 800, 500)
c8 = TCanvas('c8', 'tagged vs signal', 800, 500)
c9 = TCanvas('c9', 'tagged vs signal', 800, 500)

cleg = TCanvas('cleg', 'tagged vs signal', 400, 600)

stack1 = THStack("typeb1probeseta1", "; Probe Jet p_{T} (GeV); Events / 50(GeV)")
stack2 = THStack("typeb1probeseta2", "; Probe Jet p_{T} (GeV); Events / 50(GeV)")
stack3 = THStack("typeb1probeseta3", "; Probe Jet p_{T} (GeV); Events / 50(GeV)")

stack4 = THStack("typeb1tagseta1", "; b-tagged Jet p_{T} (GeV); Events / 50(GeV)")
stack5 = THStack("typeb1tagseta2", "; b-tagged Jet p_{T} (GeV); Events / 50(GeV)")
stack6 = THStack("typeb1tagseta3", "; b-tagged Jet p_{T} (GeV); Events / 50(GeV)")
tagrates = ROOT.TFile("plots/TBrate_Maker__PSET_"+options.cuts+".root")
ratedata = TFile(rootdir+"TBratefile"+options.set+"_PSET_"+options.cuts+".root")
ratettbar = TFile(rootdir+"TBratefilettbar_PSET_"+options.cuts+".root")

tagrateswsig = ROOT.TFile("plots/B_tagging_sigcont.root")

probeeta1data=ratedata.Get("pteta1pretag")
probeeta2data=ratedata.Get("pteta2pretag")
probeeta3data=ratedata.Get("pteta3pretag")

tageta1data=ratedata.Get("pteta1")
tageta2data=ratedata.Get("pteta2")
tageta3data=ratedata.Get("pteta3")

probeeta1mc=ratettbar.Get("pteta1pretag")
probeeta2mc=ratettbar.Get("pteta2pretag")
probeeta3mc=ratettbar.Get("pteta3pretag")

tageta1mc=ratettbar.Get("pteta1")
tageta2mc=ratettbar.Get("pteta2")
tageta3mc=ratettbar.Get("pteta3")

ptrebin = 10

probeeta1data.Rebin(ptrebin)
probeeta2data.Rebin(ptrebin)
probeeta3data.Rebin(ptrebin)

tageta1data.Rebin(ptrebin)
tageta2data.Rebin(ptrebin)
tageta3data.Rebin(ptrebin)

probeeta1mc.Rebin(ptrebin)
probeeta2mc.Rebin(ptrebin)
probeeta3mc.Rebin(ptrebin)

tageta1mc.Rebin(ptrebin)
tageta2mc.Rebin(ptrebin)
tageta3mc.Rebin(ptrebin)


probeeta1data.SetFillColor( kYellow )
probeeta2data.SetFillColor( kYellow )
probeeta3data.SetFillColor( kYellow )

tageta1data.SetFillColor( kYellow )
tageta2data.SetFillColor( kYellow )
tageta3data.SetFillColor( kYellow )

probeeta1mc.SetFillColor( kRed )
probeeta2mc.SetFillColor( kRed )
probeeta3mc.SetFillColor( kRed )

tageta1mc.SetFillColor( kRed )
tageta2mc.SetFillColor( kRed )
tageta3mc.SetFillColor( kRed )

tageta1mc=ratettbar.Get("pteta1")
tageta2mc=ratettbar.Get("pteta2")
tageta3mc=ratettbar.Get("pteta3")

treta1= tagrates.Get("tagrateeta1")
treta2= tagrates.Get("tagrateeta2")
treta3= tagrates.Get("tagrateeta3")

x = array( 'd' )
y = []
BPy = []
BPerryh = []
BPerryl = []


for eta in range(0,3):
	y.append([])
	for fittitle in fittitles:
		y[eta].append(array( 'd' ))
	BPy.append(array( 'd' ))
	BPerryh.append(array( 'd' ))
	BPerryl.append(array( 'd' ))

for j in range(0,1400):

	x.append(j)
	for eta in range(0,3):
		for ifit in range(0,len(fits)):
			y[eta][ifit].append(fits[ifit][eta].Eval(x[j]))
		BPy[eta].append(BTR[eta].Eval(x[j]))
		BPerryh[eta].append(BTR[eta].Eval(x[j])+sqrt(BTR_err[eta].Eval(x[j])))
		BPerryl[eta].append(BTR[eta].Eval(x[j])-sqrt(BTR_err[eta].Eval(x[j])))

# Create graphs of errors and ffor fittitle in fittitles:its

graphs = [] 
graphBP = []
graphBPerrh = []
graphBPerrl = []
mg = []
for eta in range(0,3):
	graphs.append([])
	for ifit in range(0,len(fits)):
		graphs[eta].append(TGraph(len(x),x,y[eta][ifit]))
		graphs[eta][ifit].SetLineColor(kBlue)
		graphs[eta][ifit].SetLineWidth(2)
	graphBP.append(TGraph(len(x),x,BPy[eta]))
	graphBP[eta].SetLineColor(kBlue)

	graphBPerrh.append(TGraph(len(x),x,BPy[eta]))
	graphBPerrl.append(TGraph(len(x),x,BPy[eta]))
	graphBPerrh[eta].SetLineColor(kBlue)
	graphBPerrl[eta].SetLineColor(kBlue)
	graphBP[eta].SetLineWidth(2)
	graphBPerrh[eta].SetLineWidth(2)
	graphBPerrl[eta].SetLineWidth(2)
	graphBPerrh[eta].SetLineStyle(2)
	graphBPerrl[eta].SetLineStyle(2)
	mg.append(TMultiGraph())
	mg[eta].Add(graphBP[eta])
	mg[eta].Add(graphBPerrh[eta])
	mg[eta].Add(graphBPerrl[eta])
#leg1.AddEntry(treta3,"Data Points","p")
leg1.AddEntry(graphBP[0],"Bifurcated polynomial fit","l")
leg1.AddEntry(graphBPerrh[0],"Fit uncertainty","l")


c1.cd()
prelim = ROOT.TLatex()
prelim.SetTextFont(42)
prelim.SetNDC()

chis = ROOT.TLatex()
chis.SetTextFont(42)
chis.SetNDC()

OFF = 1.1

SigFiles = [
ROOT.TFile(rootdir+"TBratefileweightedsignalright1300_PSET_"+options.cuts+".root"),
ROOT.TFile(rootdir+"TBratefileweightedsignalright1500_PSET_"+options.cuts+".root"),
ROOT.TFile(rootdir+"TBratefileweightedsignalright1700_PSET_"+options.cuts+".root"),
ROOT.TFile(rootdir+"TBratefileweightedsignalright1900_PSET_"+options.cuts+".root"),
ROOT.TFile(rootdir+"TBratefileweightedsignalright2100_PSET_"+options.cuts+".root"),
ROOT.TFile(rootdir+"TBratefileweightedsignalright2300_PSET_"+options.cuts+".root"),
ROOT.TFile(rootdir+"TBratefileweightedsignalright2700_PSET_"+options.cuts+".root")
]


c1.Divide(2,3)
c1.cd(1)

gPad.SetLeftMargin(0.16)
probeeta1data.Add(probeeta1mc,-1)
stack1.Add( probeeta1mc, "Hist" )
stack1.Add( probeeta1data, "Hist" )
stack1.SetMaximum(stack1.GetMaximum() * 1.2 )
stack1.Draw()
stack1.GetYaxis().SetTitleOffset(OFF)
stack1.GetXaxis().SetRangeUser(350,1200)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 8 TeV, 19.7 fb^{-1}   (0.00 < |#eta| < 0.50) }" )
c1.cd(3)
gPad.SetLeftMargin(0.16)
probeeta2data.Add(probeeta2mc,-1)
stack2.Add( probeeta2mc, "Hist" )
stack2.Add( probeeta2data, "Hist" )
stack2.SetMaximum(stack2.GetMaximum() * 1.2 )
stack2.Draw()
stack2.GetYaxis().SetTitleOffset(OFF)
stack2.GetXaxis().SetRangeUser(350,1200)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 8 TeV, 19.7 fb^{-1}   (0.50 < |#eta| < 1.15) }" )
c1.cd(5)
gPad.SetLeftMargin(0.16)
probeeta3data.Add(probeeta3mc,-1)
stack3.Add( probeeta3mc, "Hist" )
stack3.Add( probeeta3data, "Hist" )
stack3.SetMaximum(stack3.GetMaximum() * 1.2 )
stack3.Draw()
stack3.GetYaxis().SetTitleOffset(OFF)
stack3.GetXaxis().SetRangeUser(350,1200)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 8 TeV, 19.7 fb^{-1}   (1.15 < |#eta| < 2.40) }" )
c1.cd(2)
gPad.SetLeftMargin(0.16)
tageta1data.Add(tageta1mc,-1)
stack4.Add( tageta1mc, "Hist" )
stack4.Add( tageta1data, "Hist" )
stack4.SetMaximum(stack4.GetMaximum() * 1.2 )
stack4.Draw()
stack4.GetYaxis().SetTitleOffset(OFF)
stack4.GetXaxis().SetRangeUser(350,1200)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 8 TeV, 19.7 fb^{-1}   (0.00 < |#eta| < 0.50) }" )
c1.cd(4)
gPad.SetLeftMargin(0.16)
tageta2data.Add(tageta2mc,-1)
stack5.Add( tageta2mc, "Hist" )
stack5.Add( tageta2data, "Hist" )
stack5.SetMaximum(stack5.GetMaximum() * 1.2 )
stack5.Draw()
stack5.GetYaxis().SetTitleOffset(OFF)
stack5.GetXaxis().SetRangeUser(350,1200)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 8 TeV, 19.7 fb^{-1}   (0.50 < |#eta| < 1.15) }" )
c1.cd(6)
gPad.SetLeftMargin(0.16)
tageta3data.Add(tageta3mc,-1)
stack6.Add( tageta3mc, "Hist" )
stack6.Add( tageta3data, "Hist" )
stack6.SetMaximum(stack6.GetMaximum() * 1.2 )
stack6.Draw()
stack6.GetYaxis().SetTitleOffset(OFF)
stack6.GetXaxis().SetRangeUser(350,1200)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 8 TeV, 19.7 fb^{-1}   (1.15 < |#eta| < 2.40) }" )
c1.RedrawAxis()
c1.Print('plots/tagsandprobes.root', 'root')
c1.Print('plots/tagsandprobes.pdf', 'pdf')



c7.cd()
stack4.SetMinimum(0.001)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 8 TeV, 19.7 fb^{-1}   (0.00 < |#eta| < 0.50) }" )
stack4.Draw()
stack4.GetYaxis().SetTitleOffset(0.8)
c8.cd()
stack5.SetMinimum(0.001)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 8 TeV, 19.7 fb^{-1}   (0.50 < |#eta| < 1.15) }" )
stack5.Draw()
stack5.GetYaxis().SetTitleOffset(0.8)
c9.cd()
stack6.SetMinimum(0.001)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 8 TeV, 19.7 fb^{-1}   (1.15 < |#eta| < 2.40) }" )
stack6.Draw()
stack6.GetYaxis().SetTitleOffset(0.8)

leg2.AddEntry(probeeta1data,"QCD","f")
leg2.AddEntry(probeeta1mc,"ttbar","f")



for ifile in range(0,len(SigFiles)):
	if ifile<4:
		colorassn = ifile+1
	else:
		colorassn = ifile+2
		
	nseta1 = SigFiles[ifile].Get("pteta1")
	nseta2 = SigFiles[ifile].Get("pteta2")
	nseta3 = SigFiles[ifile].Get("pteta3")
	nseta1.SetLineColor(colorassn)
	nseta2.SetLineColor(colorassn)
	nseta3.SetLineColor(colorassn)
	nseta1.Rebin(ptrebin)
	nseta2.Rebin(ptrebin)
	nseta3.Rebin(ptrebin)

	scalefac = (Luminosity*Sigma[ifile]*TeffScale*Kfac)/MCentries[ifile]
	nseta1.Scale(scalefac)
	nseta2.Scale(scalefac)
	nseta3.Scale(scalefac)
	c7.cd()
	nseta1.Draw("samehist")
	c8.cd()
	nseta2.Draw("samehist")
	c9.cd()
	nseta3.Draw("samehist")
	leg2.AddEntry(nseta1,"signal at "+str(mass[ifile])+"GeV","l")
c7.SetLogy()
c8.SetLogy()
c9.SetLogy()

c7.cd()
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 8 TeV, 19.7 fb^{-1}   (0.00 < |#eta| < 0.50)  }" )
#leg2.Draw()
c8.cd()
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 8 TeV, 19.7 fb^{-1}   (0.50 < |#eta| < 1.15) }" )
#leg2.Draw()
c9.cd()
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 8 TeV, 19.7 fb^{-1}   (1.15 < |#eta| < 2.40) }" )
#leg2.Draw()

cleg.cd()
leg2.Draw()
cleg.Print('plots/legend.pdf', 'pdf')
c7.RedrawAxis()
c8.RedrawAxis()
c9.RedrawAxis()
c7.Print('plots/sigvsTReta1.root', 'root')
c7.Print('plots/sigvsTReta1.pdf', 'pdf')
c8.Print('plots/sigvsTReta2.root', 'root')
c8.Print('plots/sigvsTReta2.pdf', 'pdf')
c9.Print('plots/sigvsTReta3.root', 'root')
c9.Print('plots/sigvsTReta3.pdf', 'pdf')

trs = [
treta1,
treta2,
treta3
]

c4.cd()
c4.SetLeftMargin(0.16)

etastrings = [
'0.00 < |#eta| < 0.50',
'0.50 < |#eta| < 1.15',
'1.15 < |#eta| < 2.40'
]

for eta in range(0,3):
	for ifit in range(0,len(fits)):

		trs[eta].SetTitle(';p_{T} (GeV);Average b-tagging rate')
		trs[eta].GetYaxis().SetTitleOffset(0.8)
		trs[eta].SetMaximum(0.20)
		trs[eta].SetMinimum(0.008)
		trs[eta].SetStats(0)
		trs[eta].Draw("histe")
		graphs[eta][ifit].Draw()

		c4.RedrawAxis()
		c4.Print('plots/tagrateeta'+str(eta)+fits[ifit]+'PSET_'+options.cuts+'.root', 'root')
		c4.Print('plots/tagrateeta'+str(eta)+fits[ifit]+'PSET_'+options.cuts+'.pdf', 'pdf')
	trs[eta].SetTitle(';p_{T} (GeV);Average b-tagging rate')
	trs[eta].GetYaxis().SetTitleOffset(0.8)
	trs[eta].SetMaximum(0.20)
	trs[eta].SetMinimum(0.008)
	trs[eta].SetStats(0)
	trs[eta].Draw("histe")
	mg[eta].Draw()

	leg1.Draw()
	prelim = ROOT.TLatex()
	prelim.SetTextFont(42)
	prelim.SetNDC()

	prelim.DrawLatex( 0.2, 0.5, "#scale[1.0]{"+etastring[eta]+"}" )
	insertlogo( c4, 2, 11 )
	chis.DrawLatex( 0.20, 0.6, "#scale[1.0]{#chi^{2} / dof = "+strf(chi2eta1/ndofeta1)+"}" )
	c4.RedrawAxis()
	c4.Print('plots/tagrateeta'+str(eta)+'fitBP.root', 'root')
	c4.Print('plots/tagrateeta'+str(eta)+'fitBP.pdf', 'pdf')
	c4.Print('plots/tagrateeta'+str(eta)+'fitBP.png', 'png')




