
import os
import copy
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

parser.add_option('--batch', metavar='F', action='store_true',
                  default=False,
                  dest='batch',
                  help='batch')
(options, args) = parser.parse_args()


if options.batch:
	ROOT.gROOT.SetBatch(True)
	ROOT.PyConfig.IgnoreCommandLineOptions = True

rootdir="rootfiles/"
import Wprime_Functions	
from Wprime_Functions import *

gROOT.Macro("rootlogon.C")
#gROOT.LoadMacro("insertlogo.C+")

BTR = BTR_Init('Bifpoly',options.cuts,'',options.set)
BTR_err = BTR_Init('Bifpoly_err',options.cuts,'',options.set)

fittitles = ["pol0","pol2","pol3","FIT","expofit"]
fitlegs = ["Constant","2nd Degree Polynomial","3rd Degree Polynomial","Bifurcated Polynomial","a*(b+x)/(c+x)+d*x","Exponential + Constant"]
fits = []
for fittitle in fittitles:
	fits.append(BTR_Init(fittitle,options.cuts,'',options.set))

leg1 = TLegend(0.45,0.57,.84,.78)
leg1.SetFillColor(0)
leg1.SetBorderSize(0)

leg2 = TLegend(0.,0.,1.,1.)
leg2.SetFillColor(0)
leg2.SetBorderSize(0)

legmulti = TLegend(0.45,0.5,.75,.78)
legmulti.SetFillColor(0)
legmulti.SetBorderSize(0)

#output = ROOT.TFile( "fitting.root", "recreate" )
#output.cd()
c1 = TCanvas('c1', 'Tagrate numerator and deominator', 1000, 1300)
c4 = TCanvas('c4', '', 800, 500)

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

tagrates = ROOT.TFile("plots/TBrate_Maker_"+options.set+"_PSET_"+options.cuts+".root")

if options.set == "data":
	ratedata = TFile(rootdir+"TBratefile"+options.set+"_PSET_"+options.cuts+".root")
else:
	ratedata = TFile(rootdir+"TBratefile"+options.set+"_PSET_"+options.cuts+"weighted.root")

ratettbar = TFile(rootdir+"TBratefilettbar_PSET_"+options.cuts+"weighted.root")
ratest = TFile(rootdir+"TBratefileST_PSET_"+options.cuts+".root")
#tagrateswsig = ROOT.TFile("plots/B_tagging_sigcont.root")

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


probeeta1st=ratest.Get("pteta1pretag")
probeeta2st=ratest.Get("pteta2pretag")
probeeta3st=ratest.Get("pteta3pretag")

tageta1st=ratest.Get("pteta1")
tageta2st=ratest.Get("pteta2")
tageta3st=ratest.Get("pteta3")

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



probeeta1st.Rebin(ptrebin)
probeeta2st.Rebin(ptrebin)
probeeta3st.Rebin(ptrebin)

tageta1st.Rebin(ptrebin)
tageta2st.Rebin(ptrebin)
tageta3st.Rebin(ptrebin)



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

probeeta1st.SetFillColor( 6 )
probeeta2st.SetFillColor( 6 )
probeeta3st.SetFillColor( 6 )

tageta1st.SetFillColor( 6 )
tageta2st.SetFillColor( 6 )
tageta3st.SetFillColor( 6 )

#treta1= tagrates.Get("tagrateeta1")
#treta2= tagrates.Get("tagrateeta2")
#treta3= tagrates.Get("tagrateeta3")

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

for eta in range(0,3):
	graphs.append([])

	for ifit in range(0,len(fits)):
		graphs[eta].append(TGraph(len(x),x,y[eta][ifit]))
		graphs[eta][ifit].SetLineColor(kBlue)
		graphs[eta][ifit].SetLineWidth(2)
	graphBP.append(TGraph(len(x),x,BPy[eta]))
	graphBP[eta].SetLineColor(kBlue)

	graphBPerrh.append(TGraph(len(x),x,BPerryh[eta]))
	graphBPerrl.append(TGraph(len(x),x,BPerryl[eta]))
	graphBPerrh[eta].SetLineColor(kBlue)
	graphBPerrl[eta].SetLineColor(kBlue)
	graphBP[eta].SetLineWidth(2)
	graphBPerrh[eta].SetLineWidth(2)
	graphBPerrl[eta].SetLineWidth(2)
	graphBPerrh[eta].SetLineStyle(2)
	graphBPerrl[eta].SetLineStyle(2)




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
ROOT.TFile(rootdir+"TBratefileweightedsignalright1200_PSET_"+options.cuts+".root"),
ROOT.TFile(rootdir+"TBratefileweightedsignalright2000_PSET_"+options.cuts+".root"),
ROOT.TFile(rootdir+"TBratefileweightedsignalright2800_PSET_"+options.cuts+".root")
]


c1.Divide(2,3)
c1.cd(1)


if options.set=='data':
	probeeta3data.Add(probeeta3mc,-1)
	probeeta1data.Add(probeeta1mc,-1)
	probeeta2data.Add(probeeta2mc,-1)
	tageta1data.Add(tageta1mc,-1)
	tageta2data.Add(tageta2mc,-1)
	tageta3data.Add(tageta3mc,-1)

	probeeta3data.Add(probeeta3st,-1)
	probeeta1data.Add(probeeta1st,-1)
	probeeta2data.Add(probeeta2st,-1)
	tageta1data.Add(tageta1st,-1)
	tageta2data.Add(tageta2st,-1)
	tageta3data.Add(tageta3st,-1)


gPad.SetLeftMargin(0.16)


stack1.Add( probeeta1st, "Hist" )
stack1.Add( probeeta1mc, "Hist" )
stack1.Add( probeeta1data, "Hist" )
stack1.SetMaximum(stack1.GetMaximum() * 1.2 )
stack1.Draw()
stack1.GetYaxis().SetTitleOffset(OFF)
stack1.GetXaxis().SetRangeUser(350,1400)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 13 TeV   (0.00 < |#eta| < 0.50) }" )
c1.cd(3)
gPad.SetLeftMargin(0.16)


stack2.Add( probeeta2st, "Hist" )
stack2.Add( probeeta2mc, "Hist" )
stack2.Add( probeeta2data, "Hist" )
stack2.SetMaximum(stack2.GetMaximum() * 1.2 )
stack2.Draw()
stack2.GetYaxis().SetTitleOffset(OFF)
stack2.GetXaxis().SetRangeUser(350,1400)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 13 TeV   (0.50 < |#eta| < 1.15) }" )
c1.cd(5)
gPad.SetLeftMargin(0.16)


stack3.Add( probeeta3st, "Hist" )
stack3.Add( probeeta3mc, "Hist" )
stack3.Add( probeeta3data, "Hist" )
stack3.SetMaximum(stack3.GetMaximum() * 1.2 )
stack3.Draw()
stack3.GetYaxis().SetTitleOffset(OFF)
stack3.GetXaxis().SetRangeUser(350,1400)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 13 TeV   (1.15 < |#eta| < 2.40) }" )
c1.cd(2)
gPad.SetLeftMargin(0.16)

stack4.Add( tageta1st, "Hist" )
stack4.Add( tageta1mc, "Hist" )
stack4.Add( tageta1data, "Hist" )
stack4.SetMaximum(stack4.GetMaximum() * 1.2 )
stack4.Draw()
stack4.GetYaxis().SetTitleOffset(OFF)
stack4.GetXaxis().SetRangeUser(350,1400)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 13 TeV   (0.00 < |#eta| < 0.50) }" )
c1.cd(4)
gPad.SetLeftMargin(0.16)

stack5.Add( tageta2st, "Hist" )
stack5.Add( tageta2mc, "Hist" )
stack5.Add( tageta2data, "Hist" )
stack5.SetMaximum(stack5.GetMaximum() * 1.2 )
stack5.Draw()
stack5.GetYaxis().SetTitleOffset(OFF)
stack5.GetXaxis().SetRangeUser(350,1400)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 13 TeV   (0.50 < |#eta| < 1.15) }" )
c1.cd(6)
gPad.SetLeftMargin(0.16)

stack6.Add( tageta3st, "Hist" )
stack6.Add( tageta3mc, "Hist" )
stack6.Add( tageta3data, "Hist" )
stack6.SetMaximum(stack6.GetMaximum() * 1.2 )
stack6.Draw()
stack6.GetYaxis().SetTitleOffset(OFF)
stack6.GetXaxis().SetRangeUser(350,1400)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 13 TeV   (1.15 < |#eta| < 2.40) }" )
c1.RedrawAxis()
c1.Print('plots/tagsandprobes.root', 'root')
c1.Print('plots/tagsandprobes.pdf', 'pdf')

stacks = [stack1,stack2,stack3,stack4,stack5,stack6]
for i in range(1,7):
	stacks[i-1].SetMinimum(0.1)
	c1.cd(i).SetLogy()
c1.Print('plots/tagsandprobes_semilog.root', 'root')
c1.Print('plots/tagsandprobes_semilog.pdf', 'pdf')

c7.cd()
stack4.SetMinimum(0.001)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 13 TeV   (0.00 < |#eta| < 0.50) }" )
stack4.Draw()
stack4.GetYaxis().SetTitleOffset(0.8)
c8.cd()
stack5.SetMinimum(0.001)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 13 TeV   (0.50 < |#eta| < 1.15) }" )
stack5.Draw()
stack5.GetYaxis().SetTitleOffset(0.8)
c9.cd()
stack6.SetMinimum(0.001)
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 13 TeV   (1.15 < |#eta| < 2.40) }" )
stack6.Draw()
stack6.GetYaxis().SetTitleOffset(0.8)

leg2.AddEntry(probeeta1data,"QCD","f")
leg2.AddEntry(probeeta1mc,"ttbar","f")
leg2.AddEntry(probeeta1mc,"singletop","f")


mass = [1300,2000,2700]
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
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 13 TeV   (0.00 < |#eta| < 0.50)  }" )
#leg2.Draw()
c8.cd()
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 13 TeV   (0.50 < |#eta| < 1.15) }" )
#leg2.Draw()
c9.cd()
prelim.DrawLatex( 0.15, 0.91, "#scale[1.0]{CMS Preliminary #sqrt{s} = 13 TeV   (1.15 < |#eta| < 2.40) }" )
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




c10 = TCanvas('c10', 'QCD mass comp', 700, 500)
ModMf = ROOT.TFile("ModMassFile_PSET_"+options.cuts+".root")

ModM = ModMf.Get("bmasspost")
ModMd = ModMf.Get("bmasspre")

#ModM.Rebin(4)
#ModMd.Rebin(4)
#ModM.Scale(1/ModM.Integral())
#ModMd.Scale(1/ModMd.Integral())
ModM.SetLineColor(2)
ModMd.SetLineColor(3)
leg2.AddEntry(ModM,"post b-tagging" ,"l")
leg2.AddEntry(ModMd,"pre b-tagging" ,"l")

ModM.SetTitle(';M_{Jet} (GeV);Fraction')
ModM.SetStats(0)
gPad.SetLeftMargin(0.18)
ModM.GetYaxis().SetTitleOffset(0.95)
#ModM.SetMaximum(ModMd.GetMaximum()*1.2)
#ModM.SetMinimum(0.0)
ModM.Draw("hist")
ModMd.Draw("histsame")
leg2.Draw()
c10.Print('plots/QCDbmasscomp.pdf', 'pdf')
c10.Print('plots/QCDbmasscomp.root', 'root')
c11 = TCanvas('c11', 'Modmass plot', 700, 500)
ModM.Divide(ModM,ModMd)

#ModM.Draw("e")
ModMup = ModM.Clone()
ModMdown = ModM.Clone()
for ibin in range(0,ModM.GetNbinsX()+1):
	ModMup.SetBinContent(ibin,(1 + 0.5*(ModM.GetBinContent(ibin)-1)))
	ModMdown.SetBinContent(ibin,max(0.0,1 + 1.5*(ModM.GetBinContent(ibin)-1)))
x = []
y = []
yh = []
yl = []
x,y,yl,yh = array( 'd' ),array( 'd' ),array( 'd' ),array( 'd' )
if options.cuts=='rate_sideband3':
	start = 70.
	xarange = 200.
else:
	start = 0.
	xarange = 70.	  
for i in range(0,1000):
	x.append(start+i*xarange/1000.)
	y.append(ModM.Interpolate(x[i]))
	yh.append(ModMup.Interpolate(x[i]))	
	yl.append(ModMdown.Interpolate(x[i]))	
yg = TGraph(len(x), x, y)
yhg = TGraph(len(x), x, yh)
ylg = TGraph(len(x), x, yl)

yhg.SetLineColor(4)
ylg.SetLineColor(4)
yhg.SetLineWidth(2)
ylg.SetLineWidth(2)
yhg.SetLineStyle(2)
ylg.SetLineStyle(2)
yg.SetLineColor(1)
mmmg = ROOT.TMultiGraph()
mmmg.SetTitle(';M_{Jet} (GeV);Weight')
gPad.SetLeftMargin(0.18)

mmmg.SetMaximum(3.)
mmmg.SetMinimum(0.)
mmmg.Add(yg)
mmmg.Add(yhg)
mmmg.Add(ylg)
mmmg.Draw("al")
mmmg.GetYaxis().SetTitleOffset(0.95)
ModM.Draw("samee")
#ModMup.SetLineColor(4)
#ModMdown.SetLineColor(4)
#ModMup.SetLineWidth(2)
#ModMdown.SetLineWidth(2)
#ModMup.SetLineStyle(2)
#ModMdown.SetLineStyle(2)
#ModMup.Draw("hist same")
#ModMdown.Draw("hist same")
c11.Print('plots/ModmassPlot.pdf', 'pdf')
c11.Print('plots/ModmassPlot.root', 'root')





trs = [None]*3

trs[0]= tagrates.Get("tagrateeta1")
trs[1]= tagrates.Get("tagrateeta2")
trs[2]= tagrates.Get("tagrateeta3")

c4.cd()
c4.SetLeftMargin(0.16)

etastring = [
'0.00 < |#eta| < 0.50',
'0.50 < |#eta| < 1.15',
'1.15 < |#eta| < 2.40'
]
legf=True

mg1 = []
for eta in range(0,3):
	print 'eta '+ str(eta+1)
	trs[eta].SetTitle(';p_{T} (GeV);Average b-tagging rate')
	trs[eta].GetYaxis().SetTitleOffset(0.8)
	trs[eta].SetMaximum(0.70)
	trs[eta].SetMinimum(0.1)
	trs[eta].SetStats(0)
	c4.cd()
	trs[eta].Draw("histe")

	graphBP[eta].Draw("same")
	graphBPerrh[eta].Draw("same")
	graphBPerrl[eta].Draw("same")
	#mg[eta].Draw("same")

	#leg1.Draw()
	#prelim = ROOT.TLatex()
	#prelim.SetTextFont(42)
	#prelim.SetNDC()

	#prelim.DrawLatex( 0.2, 0.5, "#scale[1.0]{"+etastring[eta]+"}" )
	#insertlogo( c4, 2, 11 )
	#chis.DrawLatex( 0.20, 0.6, "#scale[1.0]{#chi^{2} / dof = "+strf(chi2eta1/ndofeta1)+"}" )
	leg1.Draw()
	c4.RedrawAxis()

	c4.Print('plots/tagrateeta'+str(eta+1)+'fitBP.root', 'root')
	c4.Print('plots/tagrateeta'+str(eta+1)+'fitBP.pdf', 'pdf')
	c4.Print('plots/tagrateeta'+str(eta+1)+'fitBP.png', 'png')

	c4multi = TCanvas('c4multi', '', 800, 500)
	trs[eta].Draw("histe")
	#trsmulti = trs[eta].Clone('trsmulti')

	mg1.append(TMultiGraph())

	for ifit in range(0,len(fits)):
		c4.cd()
		print ifit
		trs[eta].SetTitle(';p_{T} (GeV);Average b-tagging rate')
		trs[eta].GetYaxis().SetTitleOffset(0.8)
		trs[eta].SetMaximum(0.80)
		trs[eta].SetMinimum(0.2)
		trs[eta].SetStats(0)
		trs[eta].Draw("histe")
		graphs[eta][ifit].Draw('same')

		c4.RedrawAxis()
		c4.Print('plots/tagrateeta'+str(eta+1)+fittitles[ifit]+'PSET_'+options.cuts+'.root', 'root')
		c4.Print('plots/tagrateeta'+str(eta+1)+fittitles[ifit]+'PSET_'+options.cuts+'.pdf', 'pdf')
		
		#if ifit==0:
		#	trs[eta].Draw("histe")
		color = ifit+1
		if color>=5:
			color+=1
		graphs[eta][ifit].SetLineColor(color)
		if legf:
			legmulti.AddEntry(graphs[eta][ifit],fitlegs[ifit],'l')


		mg1[-1].Add(graphs[eta][ifit])

	legf=False
	c4multi.cd()

	mg1[-1].Draw()
	trs[eta].Draw("histesame")
	legmulti.Draw()
	c4multi.Update()
	c4multi.Print('plots/tagrateeta'+str(eta+1)+'fitmulti.root', 'root')
	c4multi.Print('plots/tagrateeta'+str(eta+1)+'fitmulti.pdf', 'pdf')
	c4multi.Print('plots/tagrateeta'+str(eta+1)+'fitmulti.png', 'png')
	c4multi.Close()
	#del mg1[-1]


tagrates.Close()
ratedata.Close()
ratettbar.Close()
#tagrateswsig.Close()
SigFiles[0].Close()
SigFiles[1].Close()
SigFiles[2].Close()





