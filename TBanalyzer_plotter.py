

import os
import array
import glob
import math
import ROOT
import sys
from ROOT import *
from array import *
from optparse import OptionParser
gROOT.Macro("rootlogon.C")
gROOT.LoadMacro("insertlogo.C+")
parser = OptionParser()

parser.add_option('-c', '--cuts', metavar='F', type='string', action='store',
                  default	=	'default',
                  dest		=	'cuts',
                  help		=	'Cuts type (ie default, rate, etc)')
parser.add_option('-s', '--set', metavar='F', type='string', action='store',
                  default	=	'data',
                  dest		=	'set',
                  help		=	'data or QCD')
parser.add_option('--batch', metavar='F', action='store_true',
                  default=False,
                  dest='batch',
                  help='batch')
(options, args) = parser.parse_args()

cuts = options.cuts

if options.batch:
	ROOT.gROOT.SetBatch(True)
	ROOT.PyConfig.IgnoreCommandLineOptions = True

text = ''
regions = ['sideband2','sideband3']
regionsstr = ['subjet b-tag inverted','b candidate softdrop mass inverted']

if options.cuts in regions:
	text = regionsstr[regions.index(options.cuts)]

		

import Wprime_Functions	
from Wprime_Functions import *

st1= ROOT.THStack("st1", "st1")

leg = TLegend(0.45, 0.4, 0.84, 0.75)
leg.SetFillColor(0)
leg.SetBorderSize(0)

rebin =4

Mult = 1.0

#Kfac=1.2



sigf = [
ROOT.TFile("rootfiles/TBanalyzerweightedsignalright1400_Trigger_nominal_none_PSET_"+options.cuts+".root"),
ROOT.TFile("rootfiles/TBanalyzerweightedsignalright2000_Trigger_nominal_none_PSET_"+options.cuts+".root"),
ROOT.TFile("rootfiles/TBanalyzerweightedsignalright2600_Trigger_nominal_none_PSET_"+options.cuts+".root")
]

sigh = [
sigf[0].Get("Mtb"),
sigf[1].Get("Mtb"),
sigf[2].Get("Mtb")
]

sigh[0].Scale(Mult)
sigh[1].Scale(Mult)
sigh[2].Scale(Mult)

sigh[0].Rebin(rebin)
sigh[1].Rebin(rebin)
sigh[2].Rebin(rebin)

sigh[0].SetLineStyle(5)
sigh[1].SetLineStyle(6)
sigh[2].SetLineStyle(7)

sigh[0].SetLineWidth(3)
sigh[1].SetLineWidth(3)
sigh[2].SetLineWidth(3)

sigh[0].SetLineColor(1)
sigh[1].SetLineColor(3)
sigh[2].SetLineColor(4)

STmc 	= ROOT.TFile("rootfiles/TBanalyzerST_Trigger_nominal_none_PSET_"+options.cuts+".root")
TTmc 	= ROOT.TFile("rootfiles/TBanalyzerttbar_none_PSET_"+options.cuts+"weighted.root")
if options.set == 'data':
	setstring = '' 
	print "running on data"
	DataB11 = ROOT.TFile("rootfiles/TBanalyzerdata_Trigger_nominal_none_PSET_"+options.cuts+".root")
	DataFS 	= DataB11.Get("Mtb") 			
	DataBE 	= DataB11.Get("QCDbkg") 		
	datapointcolor = 1	
elif options.set == 'QCD':
	setstring = 'QCD' 
	DataB11 = ROOT.TFile("rootfiles/TBanalyzerQCD_Trigger_nominal_none_PSET_"+options.cuts+".root")
	DataFS 	= DataB11.Get("Mtb") 			
	DataBE 	= DataB11.Get("QCDbkg") 		
	DataFS.Add(TTmc.Get("Mtb")) 
	datapointcolor = 4
else:
	print 'Error: Set selection invalid.'




Datamodmdown = ROOT.TFile("rootfiles/TBanalyzerdata_Trigger_nominal_none_modm_down_PSET_"+options.cuts+".root")
Datamodmup = ROOT.TFile("rootfiles/TBanalyzerdata_Trigger_nominal_none_modm_up_PSET_"+options.cuts+".root")

DataQCDmodmup = Datamodmup.Get("QCDbkg")
DataQCDmodmdown = Datamodmdown.Get("QCDbkg")
DataBE2d = DataB11.Get("QCDbkg2D") 

c1 = TCanvas('c1', 'QCD Full selection vs b pt tagging background', 700, 600)



main = ROOT.TPad("main", "main", 0, 0.3, 1, 1)
sub = ROOT.TPad("sub", "sub", 0, 0, 1, 0.3)

main.SetLeftMargin(0.16)
main.SetRightMargin(0.05)
main.SetTopMargin(0.11)
main.SetBottomMargin(0.0)

sub.SetLeftMargin(0.16)
sub.SetRightMargin(0.05)
sub.SetTopMargin(0)
sub.SetBottomMargin(0.3)

main.Draw()
sub.Draw()

main.cd()


TTmcPtScaleUp = ROOT.TFile("rootfiles/TBanalyzerttbar_JES_up_PSET_"+options.cuts+"weighted.root")
TTmcPtScaleDown = ROOT.TFile("rootfiles/TBanalyzerttbar_JES_down_PSET_"+options.cuts+"weighted.root")


TTmcPtSmearUp = ROOT.TFile("rootfiles/TBanalyzerttbar_JER_up_PSET_"+options.cuts+"weighted.root")
TTmcPtSmearDown = ROOT.TFile("rootfiles/TBanalyzerttbar_JER_down_PSET_"+options.cuts+"weighted.root")


TTmcQ2ScaleUp = ROOT.TFile("rootfiles/TBanalyzerttbarscaleup_Trigger_nominal_none_PSET_"+options.cuts+"weighted.root")
TTmcQ2ScaleDown = ROOT.TFile("rootfiles/TBanalyzerttbarscaledown_Trigger_nominal_none_PSET_"+options.cuts+"weighted.root")




#TTmcEtaSmearUp = ROOT.TFile("rootfiles/TBanalyzerttbar_Trigger_nominal_EtaSmearUp_PSET_"+options.cuts+".root")
#TTmcEtaSmearDown = ROOT.TFile("rootfiles/TBanalyzerttbar_Trigger_nominal_EtaSmearDown_PSET_"+options.cuts+".root")

output = ROOT.TFile( "TBanalyzer_output_PSET_"+options.cuts+".root", "recreate" )
output.cd()

TTmcFS = TTmc.Get("Mtb")
TTmcBE = TTmc.Get("QCDbkg")
TTmcBE2d = TTmc.Get("QCDbkg2D")
TTmcBEh = TTmc.Get("QCDbkgh")
TTmcBEl = TTmc.Get("QCDbkgl")

STmcFS = STmc.Get("Mtb")
STmcBE = STmc.Get("QCDbkg")
STmcBE2d = STmc.Get("QCDbkg2D")
STmcBEh = STmc.Get("QCDbkgh")
STmcBEl = STmc.Get("QCDbkgl")

DataBEh = DataB11.Get("QCDbkgh")
DataBEl = DataB11.Get("QCDbkgl")

TTmcFSScaleUp = TTmcPtScaleUp.Get("Mtb")
TTmcFSScaleDown = TTmcPtScaleDown.Get("Mtb")

TTmcFSQ2ScaleUp = TTmcQ2ScaleUp.Get("Mtb")
TTmcFSQ2ScaleDown = TTmcQ2ScaleDown.Get("Mtb")

TTmcFSPtSmearUp = TTmcPtSmearUp.Get("Mtb")
TTmcFSPtSmearDown = TTmcPtSmearDown.Get("Mtb")

#TTmcFSEtaSmearUp = TTmcEtaSmearUp.Get("Mtb")
#TTmcFSEtaSmearDown = TTmcEtaSmearDown.Get("Mtb")

TTmcFSTriggerUp = TTmc.Get("Mtbtrigup")
TTmcFSTriggerDown = TTmc.Get("Mtbtrigdown")

TTmcFSBUp =  TTmc.Get("MtbBup")
TTmcFSBDown =  TTmc.Get("MtbBdown")

TTmcFSTUp =  TTmc.Get("MtbTup")
TTmcFSTDown =  TTmc.Get("MtbTdown")

TTmcFSBUp.Rebin(rebin)
TTmcFSBDown.Rebin(rebin)

TTmcFSTUp.Rebin(rebin)
TTmcFSTDown.Rebin(rebin)

#TTmcFSptup.Rebin(rebin)
#TTmcFSptdown.Rebin(rebin)

STmcFSBUp =  STmc.Get("MtbBup")
STmcFSBDown =  STmc.Get("MtbBdown")

STmcFSBUp.Rebin(rebin)
STmcFSBDown.Rebin(rebin)


TTmcFSQ2ScaleUp.Rebin(rebin)
TTmcFSQ2ScaleDown.Rebin(rebin)

TTmcFSScaleUp.Rebin(rebin)
TTmcFSScaleDown.Rebin(rebin)

TTmcFSPtSmearUp.Rebin(rebin)
TTmcFSPtSmearDown.Rebin(rebin)

#TTmcFSEtaSmearUp.Rebin(rebin)
#TTmcFSEtaSmearDown.Rebin(rebin)

TTmcFSTriggerUp.Rebin(rebin)
TTmcFSTriggerDown.Rebin(rebin)


DataBE2d.Rebin(rebin)
TTmcFS.Rebin(rebin)
STmcFS.Rebin(rebin)
DataBE.Rebin(rebin)
DataQCDmodmup.Rebin(rebin)
DataQCDmodmdown.Rebin(rebin)

print DataFS.Integral()
DataFS.Rebin(rebin)

DataBEl.Rebin(rebin)
DataBEh.Rebin(rebin)

TTmcBE.Rebin(rebin)
TTmcBE2d.Rebin(rebin)
TTmcBEh.Rebin(rebin)
TTmcBEl.Rebin(rebin)

STmcBE.Rebin(rebin)
STmcBE2d.Rebin(rebin)
STmcBEh.Rebin(rebin)
STmcBEl.Rebin(rebin)

#subtract weighted TT pretags

unsubbkg = DataBE.Clone()

# Uncomment when using data -JL #
if options.set=='data':
	DataBE.Add(TTmcBE,-1)
	DataBE2d.Add(TTmcBE2d,-1)
	DataBEl.Add(TTmcBE,-1)
	DataBEh.Add(TTmcBE,-1)
	DataQCDmodmup.Add(TTmcBE,-1)
	DataQCDmodmdown.Add(TTmcBE,-1)

	DataBE.Add(STmcBE,-1)
	DataBE2d.Add(STmcBE2d,-1)
	DataBEl.Add(STmcBE,-1)
	DataBEh.Add(STmcBE,-1)
	DataQCDmodmup.Add(STmcBE,-1)
	DataQCDmodmdown.Add(STmcBE,-1)

output.cd()

fittitles = ["pol0","pol2","pol3","Bifpoly","expofit"]
QCDbkg_ARR = []
for ihist in range(0,len(fittitles)):
	QCDbkg_ARR.append(DataB11.Get("QCDbkg"+str(fittitles[ihist])).Rebin(rebin))

BEfiterrh = Fit_Uncertainty(QCDbkg_ARR)

output.cd()
DataQCDBEH=DataBE.Clone("DataQCDBEH")
DataQCDBEL=DataBE.Clone("DataQCDBEL")
DataTOTALBEH=DataBE.Clone("DataTOTALBEH")
DataTOTALBEL=DataBE.Clone("DataTOTALBEL")


for ibin in range(0,DataBE.GetNbinsX()):
	PtScaleup=(TTmcFSScaleUp.GetBinContent(ibin) -TTmcFS.GetBinContent(ibin))
	Q2Scaleup=(TTmcFSQ2ScaleUp.GetBinContent(ibin) -TTmcFS.GetBinContent(ibin))
	PtSmearup=(TTmcFSPtSmearUp.GetBinContent(ibin) -TTmcFS.GetBinContent(ibin))
	#EtaSmearup=(TTmcFSEtaSmearUp.GetBinContent(ibin) -TTmcFS.GetBinContent(ibin))
	Triggerup=(TTmcFSTriggerUp.GetBinContent(ibin) -TTmcFS.GetBinContent(ibin))
	Btaggingup=(TTmcFSBUp.GetBinContent(ibin)-TTmcFS.GetBinContent(ibin))
	Ttaggingup=(TTmcFSTUp.GetBinContent(ibin)-TTmcFS.GetBinContent(ibin))
	#ptup=(TTmcFSptup.GetBinContent(ibin)-TTmcFS.GetBinContent(ibin))

	PtScaledown=(TTmcFSScaleDown.GetBinContent(ibin) -TTmcFS.GetBinContent(ibin))
	Q2Scaledown=(TTmcFSQ2ScaleDown.GetBinContent(ibin) -TTmcFS.GetBinContent(ibin))
	PtSmeardown=(TTmcFSPtSmearDown.GetBinContent(ibin) -TTmcFS.GetBinContent(ibin))
	#EtaSmeardown=(TTmcFSEtaSmearDown.GetBinContent(ibin) -TTmcFS.GetBinContent(ibin))
	Triggerdown=(TTmcFSTriggerDown.GetBinContent(ibin) -TTmcFS.GetBinContent(ibin))
	Btaggingdown=(TTmcFSBDown.GetBinContent(ibin)-TTmcFS.GetBinContent(ibin))
	Ttaggingdown=(TTmcFSTDown.GetBinContent(ibin)-TTmcFS.GetBinContent(ibin))
	#ptdown=(TTmcFSptdown.GetBinContent(ibin)-TTmcFS.GetBinContent(ibin))
	tsigup = 0.048*TTmcFS.GetBinContent(ibin)
	tsigdown = -0.055*TTmcFS.GetBinContent(ibin)
	ttsf = 0.215*TTmcFS.GetBinContent(ibin)
	tlumi = 0.027*TTmcFS.GetBinContent(ibin)

	ups = [PtScaleup,Q2Scaleup,PtSmearup,Triggerup,Btaggingup,ttsf,tsigup,tlumi]
	downs = [PtScaledown,Q2Scaledown,PtSmeardown,Triggerdown,Btaggingdown,-1*ttsf,tsigdown,-1*tlumi]
	
	sigsqup = 0.
	sigsqdown = 0.

	for i in range(0,len(ups)):
		upsig = max(ups[i],downs[i],0.)
		downsig = min(ups[i],downs[i],0.)
		sigsqup+=upsig*upsig
		sigsqdown+=downsig*downsig

	#CrossSection=0.19*TTmcFS.GetBinContent(ibin)
	TTstat=TTmcFS.GetBinError(ibin)
	STstat=STmcFS.GetBinError(ibin)
	if DataBE.GetBinContent(ibin)>0:
		QCDstat=DataBE.GetBinError(ibin)
	else:
		QCDstat=0.
	QCDfit=abs(BEfiterrh.GetBinContent(ibin))
	QCDfit1=abs((DataBEh.GetBinContent(ibin)-DataBEl.GetBinContent(ibin))/2.0)
	QCDfit2=abs(DataBE2d.GetBinContent(ibin)-DataBE.GetBinContent(ibin))
	QCDmodm=abs((DataQCDmodmup.GetBinContent(ibin)-DataQCDmodmdown.GetBinContent(ibin))/2.0)
	print ""
	print ibin
	print QCDfit
	print QCDfit1
	print QCDfit2

	QCDsys=math.sqrt(QCDfit*QCDfit + QCDfit1*QCDfit1 + QCDfit2*QCDfit2 +QCDmodm*QCDmodm)
	QCDerror=math.sqrt(QCDstat*QCDstat+QCDsys*QCDsys)


	TTerrorup=math.sqrt(sigsqup+TTstat*TTstat+STstat*STstat)
	TTerrordown=math.sqrt(sigsqdown+TTstat*TTstat+STstat*STstat)
	Totalerrorup=math.sqrt(QCDerror*QCDerror+TTerrorup*TTerrorup)
	Totalerrordown=math.sqrt(QCDerror*QCDerror+TTerrordown*TTerrordown)
	DataQCDBEH.SetBinContent(ibin,DataQCDBEH.GetBinContent(ibin)+QCDerror)
	DataQCDBEL.SetBinContent(ibin,DataQCDBEL.GetBinContent(ibin)-QCDerror)
	DataTOTALBEH.SetBinContent(ibin,DataTOTALBEH.GetBinContent(ibin)+Totalerrorup)
	DataTOTALBEL.SetBinContent(ibin,DataTOTALBEL.GetBinContent(ibin)-Totalerrordown)

print "QCD total error"
print (DataQCDBEH.Integral()-DataBE.Integral())/DataBE.Integral()
print 
DataQCDBEH.Write()
DataQCDBEL.Write()

DataTOTALBEH.Write()
DataTOTALBEL.Write()

DataTOTALBEL.Add(TTmcFS)
DataTOTALBEH.Add(TTmcFS)

DataTOTALBEL.Add(STmcFS)
DataTOTALBEH.Add(STmcFS)		

DataBE.SetFillColor(kYellow)
TTmcFS.SetFillColor(kRed)
STmcFS.SetFillColor(6)

DataTOTALBEH.SetLineColor(kBlue)
DataTOTALBEH.SetLineWidth(2)
#DataTOTALBEH.SetLineStyle(2)

centerqcd = DataTOTALBEL.Clone("centerqcd")
centerqcd.SetFillColor(kYellow)
centerqcd.Add(TTmcFS,-1)
centerqcd.Add(STmcFS,-1)

DataTOTALBEL.SetLineColor(kBlue)
DataTOTALBEL.SetLineWidth(2)
#DataTOTALBEL.SetFillColor(0)
#DataTOTALBEL.SetLineStyle(0)
#DataTOTALBEL.SetLineWidth(2)
#DataTOTALBEL.SetLineStyle(2)

sigst= ROOT.THStack("sigst", "sigst")
sigma = DataTOTALBEH.Clone("sigma")
sigma.SetFillStyle(3245)
sigma.SetFillColor(1)
sigma.SetLineColor(0)
centerqcd.SetLineColor(kYellow)

sigma.Add(DataTOTALBEL,-1)
sigst.Add(STmcFS)
sigst.Add(TTmcFS)
sigst.Add(centerqcd)
sigst.Add(sigma)

st1.Add(STmcFS)
st1.Add(TTmcFS)
st1.Add(DataBE)

bkgline=st1.GetStack().Last().Clone("bkgline")
bkgline.SetFillColor(0)
bkgline.SetFillStyle(0)
DataFS.SetLineColor(datapointcolor)
DataFS.SetMarkerColor(datapointcolor)
#leg.AddEntry( DataFS, 'Data', 'P')
if options.set == 'QCD':
	leg.AddEntry( DataFS, 'QCD FS + t#bar{t} + Single top prediction', 'P')
	leg.AddEntry( DataBE, 'QCD Background', 'F')
elif options.set == 'data':
	leg.AddEntry( DataFS, 'data', 'P')
	leg.AddEntry( DataBE, 'QCD', 'F')
leg.AddEntry( TTmcFS, 't#bar{t}', 'F')
leg.AddEntry( STmcFS, 'Single top', 'F')

leg.AddEntry( sigma, '1 #sigma background uncertainty', 'F')
leg.AddEntry( sigh[0], 'W`_{R} at 1400 GeV', 'L')
leg.AddEntry( sigh[1], 'W`_{R} at 2000 GeV', 'L')
leg.AddEntry( sigh[2], 'W`_{R} at 2600 GeV', 'L')

#c1.cd()
#c1.SetLeftMargin(0.17)
#st1.GetXaxis().SetRangeUser(0,3000)
st1.SetMaximum(DataTOTALBEH.GetMaximum() * 1.3)
st1.SetMinimum(1.0)
st1.SetTitle(';M_{tb} (GeV);Counts per 100 GeV')
st1.Draw("hist")
gPad.SetLeftMargin(.16)
st1.GetYaxis().SetTitleOffset(0.9)
#DataTOTALBEH.Draw("histsame")
#DataTOTALBEL.Draw("histsame")
sigst.Draw("samehist")
bkgline.Draw("samehist")
sigh[0].Draw("samehist")
sigh[1].Draw("samehist")
sigh[2].Draw("samehist")


DataFS1	    = TH1D("DataFS1",     "mass W' in b+1",     	  	      140, 500, 4000 )
DataFS1.Rebin(rebin)
for ibin in range(1,DataFS.GetNbinsX()+1):
	DataFS1.SetBinContent(ibin,DataFS.GetBinContent(ibin))

DataFS1.SetBinErrorOption(DataFS1.kPoisson)
DataFS1 = DataFS
DataFS1.Draw("samepE")

leg.Draw()
prelim = TLatex()
prelim.SetNDC()


#insertlogo( main, 2, 11 )


prelim.DrawLatex( 0.5, 0.91, "#scale[0.8]{CMS Preliminary, 13 TeV, 2553 pb^{-1}}" )
prelim.DrawLatex( 0.2, 0.83, "#scale[0.8]{"+text+"}" )
sub.cd()
gPad.SetLeftMargin(.16)
totalH = st1.GetStack().Last().Clone("totalH")
#totalH.Add(TTmcFS)
pull = Make_Pull_plot( DataFS1,totalH,DataTOTALBEH,DataTOTALBEL )


	
#pull.GetXaxis().SetRangeUser(0,3000)
pull.SetFillColor(kBlue)
pull.SetTitle(';M_{tb} (GeV);(Data-Bkg)/#sigma')
pull.SetStats(0)


pull.GetYaxis().SetRangeUser(-2.9,2.9)
pull.GetXaxis().SetLabelSize(0.05)
pull.GetYaxis().SetLabelSize(0.05)


LS = .13

pull.GetYaxis().SetTitleOffset(0.4)
pull.GetXaxis().SetTitleOffset(0.9)
pull.SetStats(0)
    

pull.GetYaxis().SetLabelSize(LS)
pull.GetYaxis().SetTitleSize(LS)
pull.GetYaxis().SetNdivisions(306)
pull.GetXaxis().SetLabelSize(LS)
pull.GetXaxis().SetTitleSize(LS)

pull.Draw("hist")

line2=ROOT.TLine(500.0,0.0,4000.0,0.0)
line2.SetLineColor(0)
line1=ROOT.TLine(500.0,0.0,4000.0,0.0)
line1.SetLineStyle(2)

line2.Draw()
line1.Draw()
gPad.Update()

main.RedrawAxis()

c1.Print('plots/MtbvsBkg_BifPoly_fit_'+setstring+'PSET_'+options.cuts+'.root', 'root')
c1.Print('plots/MtbvsBkg_BifPoly_fit_'+setstring+'PSET_'+options.cuts+'.pdf', 'pdf')
c1.Print('plots/MtbvsBkg_BifPoly_fit_'+setstring+'PSET_'+options.cuts+'.png', 'png')
main.SetLogy()
st1.SetMaximum( DataBEh.GetMaximum() * 5000 )
st1.SetMinimum( 0.5)
main.RedrawAxis()

c1.Print('plots/MtbvsBkgsemilog_BifPoly_fit_'+setstring+'PSET_'+options.cuts+'.root', 'root')
c1.Print('plots/MtbvsBkgsemilog_BifPoly_fit_'+setstring+'PSET_'+options.cuts+'.pdf', 'pdf')
c1.Print('plots/MtbvsBkgsemilog_BifPoly_fit_'+setstring+'PSET_'+options.cuts+'.png', 'png')
datatotal = DataFS1.Clone("datatotal")
bkgtotal = st1.GetStack().Last().Clone("bkgtotal")
datatotal.GetXaxis().SetRangeUser(800.,2000.)
bkgtotal.GetXaxis().SetRangeUser(800.,2000.)
if options.set=='data':
	print datatotal.Chi2Test(bkgtotal,"UW")
else:
	print datatotal.Chi2Test(bkgtotal,"WW")




