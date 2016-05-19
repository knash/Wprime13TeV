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

parser.add_option('-c', '--cuts', metavar='F', type='string', action='store',
                  default	=	'default',
                  dest		=	'cuts',
                  help		=	'Cuts type (ie default, rate, etc)')
(options, args) = parser.parse_args()


def Zero(hist):
	for ibin in range(0,hist.GetXaxis().GetNbins()+1):
		hist.SetBinContent(ibin,max(0.0,hist.GetBinContent(ibin)))
	
cuts = options.cuts
gROOT.Macro("rootlogon.C")
import Wprime_Functions	
from Wprime_Functions import *

#titles=["PtScaling","TriggerWeighting","PtSmearing"]
#modsup=["ScaleUp","none","PtSmearUp"]
#modsdown=["ScaleDown","none","PtSmearDown"]
#trigsup=["nominal","up","nominal"]
#trigsdown=["nominal","down","nominal"]


LabelsU=['__jes__','__trig__','__ptsmear__']
mass = [1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900]
rebin =2
for coup in ['right']:

	output = ROOT.TFile( "limitsetting/theta/Limits_allhadronic_"+coup+"_PSET_"+cuts+".root", "recreate" )
	output.cd()


	Data = ROOT.TFile("rootfiles/TBanalyzerdata_Trigger_nominal_none_PSET_"+options.cuts+".root")
	Datamodmdown = ROOT.TFile("rootfiles/TBanalyzerdata_Trigger_nominal_none_modm_down_PSET_"+options.cuts+".root")
	Datamodmup = ROOT.TFile("rootfiles/TBanalyzerdata_Trigger_nominal_none_modm_up_PSET_"+options.cuts+".root")


	TTmc 	= ROOT.TFile("rootfiles/TBanalyzerttbar_none_PSET_"+options.cuts+"weighted.root")
	STmc = ROOT.TFile("rootfiles/TBanalyzerST_Trigger_nominal_none_PSET_"+options.cuts+".root")

	TTmcPtScaleUp = ROOT.TFile("rootfiles/TBanalyzerttbar_JES_up_PSET_"+options.cuts+"weighted.root")
	TTmcPtScaleDown = ROOT.TFile("rootfiles/TBanalyzerttbar_JES_down_PSET_"+options.cuts+"weighted.root")


	TTmcPtSmearUp = ROOT.TFile("rootfiles/TBanalyzerttbar_JER_up_PSET_"+options.cuts+"weighted.root")
	TTmcPtSmearDown = ROOT.TFile("rootfiles/TBanalyzerttbar_JER_down_PSET_"+options.cuts+"weighted.root")


	TTmcQ2ScaleUp = ROOT.TFile("rootfiles/TBanalyzerttbarscaleup_Trigger_nominal_none_PSET_"+options.cuts+"weighted.root")
	TTmcQ2ScaleDown = ROOT.TFile("rootfiles/TBanalyzerttbarscaledown_Trigger_nominal_none_PSET_"+options.cuts+"weighted.root")

	TTmcPileUp = ROOT.TFile("rootfiles/TBanalyzerttbar__pileup_up_PSET_"+options.cuts+"weighted.root")
	TTmcPileDown = ROOT.TFile("rootfiles/TBanalyzerttbar__pileup_down_PSET_"+options.cuts+"weighted.root")

	TTmcPDFUp = ROOT.TFile("rootfiles/TBanalyzerttbar_none_pdf__up_PSET_"+options.cuts+"weighted.root")
	TTmcPDFDown = ROOT.TFile("rootfiles/TBanalyzerttbar_none_pdf__down_PSET_"+options.cuts+"weighted.root")



	#STsmc = ROOT.TFile("rootfiles/TBanalyzerweightedsingletop_s_Trigger_nominal_none_PSET_"+cuts+".root")


	#STsmc = ROOT.TFile("rootfiles/TBanalyzerweightedsingletop_s_Trigger_nominal_none_PSET_"+cuts+".root")
	#STsBmc = ROOT.TFile("rootfiles/TBanalyzerweightedsingletop_sB_Trigger_nominal_none_PSET_"+cuts+".root")
	#STtmc = ROOT.TFile("rootfiles/TBanalyzerweightedsingletop_t_Trigger_nominal_none_PSET_"+cuts+".root")
	#STtBmc = ROOT.TFile("rootfiles/TBanalyzerweightedsingletop_tB_Trigger_nominal_none_PSET_"+cuts+".root")
	#STtWmc = ROOT.TFile("rootfiles/TBanalyzerweightedsingletop_tW_Trigger_nominal_none_PSET_"+cuts+".root")
	#STtWBmc = ROOT.TFile("rootfiles/TBanalyzerweightedsingletop_tWB_Trigger_nominal_none_PSET_"+cuts+".root")

	#starray = [
	#STsmc,
	#STsBmc, 
	#STtmc,
	#STtBmc,
	#STtWmc,
	#STtWBmc 
	#]

	#STsmcFS = STsmc.Get("Mtb")
	#STsBmcFS = STsBmc.Get("Mtb")
	#STtmcFS = STtmc.Get("Mtb")
	#STtBmcFS = STtBmc.Get("Mtb")
	#STtWmcFS = STtWmc.Get("Mtb")
	#STtWBmcFS = STtWBmc.Get("Mtb")

	#STsmcFS.Add(STsBmcFS)
	#STtmcFS.Add(STtBmcFS) 
	#STtWmcFS.Add(STtWBmcFS) 


	TTmcFS = TTmc.Get("Mtb")
	STmcFS = STmc.Get("Mtb")

	TTmcQCD = TTmc.Get("QCDbkg")
	TTmcQCD2d = TTmc.Get("QCDbkg2D")

	STmcQCD = STmc.Get("QCDbkg")
	STmcQCD2d = STmc.Get("QCDbkg2D")


	TTmcFSPtScaleUp = TTmcPtScaleUp.Get("Mtb")
	TTmcFSPtScaleDown = TTmcPtScaleDown.Get("Mtb")

	TTmcFSQ2ScaleUp = TTmcQ2ScaleUp.Get("Mtb")
	TTmcFSQ2ScaleDown = TTmcQ2ScaleDown.Get("Mtb")

	TTmcFSPtSmearUp = TTmcPtSmearUp.Get("Mtb")
	TTmcFSPtSmearDown = TTmcPtSmearDown.Get("Mtb")

	TTmcFSPileUp = TTmcPileUp.Get("Mtb")
	TTmcFSPileDown = TTmcPileDown.Get("Mtb")

	TTmcFSPDFUp = TTmcPDFUp.Get("Mtb")
	TTmcFSPDFDown = TTmcPDFDown.Get("Mtb")

	TTmcFSTriggerUp = TTmc.Get("Mtbtrigup")
	TTmcFSTriggerDown = TTmc.Get("Mtbtrigdown")

	TTmcFSBup =  TTmc.Get("MtbBup")
	TTmcFSBdown =  TTmc.Get("MtbBdown")

	TTmcFSTup =  TTmc.Get("MtbTup")
	TTmcFSTdown =  TTmc.Get("MtbTdown")

	DataFS = Data.Get("Mtb")
	DataQCD = Data.Get("QCDbkg")
	DataQCD2d = Data.Get("QCDbkg2D")
	DataQCDUp = Data.Get("QCDbkgh")
	DataQCDDown = Data.Get("QCDbkgl")
	DataQCDmodmup = Datamodmup.Get("QCDbkg")
	DataQCDmodmdown = Datamodmdown.Get("QCDbkg")
	#TTmcQCD2d = TTmc.Get("QCDbkg2D")
	#DataFS.Add(TTmc.Get("Mtb"))

	DataQCD.Add(TTmcQCD,-1)
	DataQCDUp.Add(TTmcQCD,-1)
	DataQCDDown.Add(TTmcQCD,-1)
	DataQCD2d.Add(TTmcQCD2d,-1)
	DataQCDmodmup.Add(TTmcQCD,-1)
	DataQCDmodmdown.Add(TTmcQCD,-1)


	DataQCD.Add(STmcQCD,-1)
	DataQCDUp.Add(STmcQCD,-1)
	DataQCDDown.Add(STmcQCD,-1)
	DataQCD2d.Add(STmcQCD2d,-1)
	DataQCDmodmup.Add(STmcQCD,-1)
	DataQCDmodmdown.Add(STmcQCD,-1)

	Zero(DataQCD)
	Zero(DataQCDUp)
	Zero(DataQCDDown)
	Zero(DataQCD2d)
	Zero(DataQCDmodmup)
	Zero(DataQCDmodmdown)
	#for f in starray:
	#	DataQCD.Add(f.Get("QCDbkg"),-1)
	#	DataQCDUp.Add(f.Get("QCDbkg"),-1)
	#	DataQCDDown.Add(f.Get("QCDbkg"),-1)
	#	DataQCD2d.Add(f.Get("QCDbkg2D"),-1)


	fittitles = ["pol0","pol2","pol3","FIT","Bifpoly","expofit"]
	QCDbkg_ARR = []
	for ihist in range(0,len(fittitles)):
		QCDbkg_ARR.append(Data.Get("QCDbkg"+str(fittitles[ihist])).Rebin(rebin))

	BEfiterrh = Fit_Uncertainty(QCDbkg_ARR)


	DataQCD2d.Rebin(rebin)
	DataFS.Rebin(rebin)
	DataQCD.Rebin(rebin)
	DataQCDmodmup.Rebin(rebin)
	DataQCDmodmdown.Rebin(rebin)


	DataQCDUp.Rebin(rebin)
	DataQCDDown.Rebin(rebin)
	#DataBEtptl.Rebin(rebin)
	#DataBEtpth.Rebin(rebin)
	TTmcFS.Rebin(rebin)
	STmcFS.Rebin(rebin)



	TTmcFSQ2ScaleUp.Rebin(rebin)
	TTmcFSQ2ScaleDown.Rebin(rebin)

	TTmcFSPtScaleUp.Rebin(rebin)
	TTmcFSPtScaleDown.Rebin(rebin)

	TTmcFSPtSmearUp.Rebin(rebin)
	TTmcFSPtSmearDown.Rebin(rebin)

	TTmcFSPileUp.Rebin(rebin)
	TTmcFSPileDown.Rebin(rebin)

	TTmcFSPDFUp.Rebin(rebin)
	TTmcFSPDFDown.Rebin(rebin)

	TTmcFSTriggerUp.Rebin(rebin)
	TTmcFSTriggerDown.Rebin(rebin)


	TTmcFSBup.Rebin(rebin)
	TTmcFSBdown.Rebin(rebin)

	TTmcFSTup.Rebin(rebin)
	TTmcFSTdown.Rebin(rebin)


	#STsmcFS.Rebin(rebin)
	#STtmcFS.Rebin(rebin)
	#STtWmcFS.Rebin(rebin)

	DataQCDE1Up = DataQCD.Clone()	
	DataQCDE2Up = DataQCD.Clone()	
	DataQCDE1Down = DataQCD.Clone()	
	DataQCDE2Down = DataQCD.Clone()
	for ibin in range(0,DataQCD.GetNbinsX()+1):
		QCDfit3=abs(DataQCD2d.GetBinContent(ibin)-DataQCD.GetBinContent(ibin))
		QCDfit2=abs(BEfiterrh.GetBinContent(ibin))
		DataQCDE1Up.SetBinContent(ibin,max(0.0,DataQCD.GetBinContent(ibin)+QCDfit3))
		DataQCDE1Down.SetBinContent(ibin,max(0.0,DataQCD.GetBinContent(ibin)-QCDfit3))
		DataQCDE2Up.SetBinContent(ibin,max(0.0,DataQCD.GetBinContent(ibin)+QCDfit2))
		DataQCDE2Down.SetBinContent(ibin,max(0.0,DataQCD.GetBinContent(ibin)-QCDfit2))
	



	DataFS.SetName("mtb_allhad__DATA")
	DataQCD.SetName("mtb_allhad__qcd")
	DataQCDUp.SetName("mtb_allhad__qcd__Fit__plus")
	DataQCDDown.SetName("mtb_allhad__qcd__Fit__minus")
	DataQCDmodmup.SetName("mtb_allhad__qcd__modm__plus")
	DataQCDmodmdown.SetName("mtb_allhad__qcd__modm__minus")
	DataQCDE1Up.SetName("mtb_allhad__qcd__TwoD__plus")
	DataQCDE1Down.SetName("mtb_allhad__qcd__TwoD__minus")
	DataQCDE2Up.SetName("mtb_allhad__qcd__Alt__plus")
	DataQCDE2Down.SetName("mtb_allhad__qcd__Alt__minus")

	TTmcFS.SetName("mtb_allhad__ttbar")
	STmcFS.SetName("mtb_allhad__st")

	TTmcFSPtScaleUp.SetName("mtb_allhad__ttbar__jes__plus")
	TTmcFSPtScaleDown.SetName("mtb_allhad__ttbar__jes__minus")

	TTmcFSQ2ScaleUp.SetName("mtb_allhad__ttbar__q2__plus")
	TTmcFSQ2ScaleDown.SetName("mtb_allhad__ttbar__q2__minus")

	TTmcFSBup.SetName("mtb_allhad__ttbar__btag__plus")
	TTmcFSBdown.SetName("mtb_allhad__ttbar__btag__minus")

	TTmcFSTup.SetName("mtb_allhad__ttbar__ttag__plus")
	TTmcFSTdown.SetName("mtb_allhad__ttbar__ttag__minus")


	#TTmcptup.SetName("mtb_allhad__ttbar__pt__plus")
	#TTmcptdown.SetName("mtb_allhad__ttbar__pt__minus")

	TTmcFSPtSmearUp.SetName("mtb_allhad__ttbar__jer__plus")
	TTmcFSPtSmearDown.SetName("mtb_allhad__ttbar__jer__minus")



	TTmcFSPileUp.SetName("mtb_allhad__ttbar__pile__plus")
	TTmcFSPileDown.SetName("mtb_allhad__ttbar__pile__minus")


	TTmcFSPDFUp.SetName("mtb_allhad__ttbar__pdf__plus")
	TTmcFSPDFDown.SetName("mtb_allhad__ttbar__pdf__minus")


	TTmcFSTriggerUp.SetName("mtb_allhad__ttbar__trig__plus")
	TTmcFSTriggerDown.SetName("mtb_allhad__ttbar__trig__minus")

	#STsmcFS.SetName("mtb_allhad__sts")
	#STtmcFS.SetName("mtb_allhad__stt")
	#STtWmcFS.SetName("mtb_allhad__sttw")

	DataFS.SetTitle("mtb_allhad__DATA")
	DataQCD.SetTitle("mtb_allhad__qcd")
	DataQCDUp.SetTitle("mtb_allhad__qcd__Fit__plus")
	DataQCDDown.SetTitle("mtb_allhad__qcd__Fit__minus")
	DataQCDE1Up.SetTitle("mtb_allhad__qcd__TwoD__plus")
	DataQCDE1Down.SetTitle("mtb_allhad__qcd__TwoD__minus")
	DataQCDE2Up.SetTitle("mtb_allhad__qcd__Alt__plus")
	DataQCDE2Down.SetTitle("mtb_allhad__qcd__Alt__minus")
	DataQCDmodmup.SetTitle("mtb_allhad__qcd__modm__plus")
	DataQCDmodmdown.SetTitle("mtb_allhad__qcd__modm__minus")



	TTmcFS.SetTitle("mtb_allhad__ttbar")
	STmcFS.SetTitle("mtb_allhad__st")
	TTmcFSPtScaleUp.SetTitle("mtb_allhad__ttbar__jes__plus")
	TTmcFSPtScaleDown.SetTitle("mtb_allhad__ttbar__jes__minus")
	
	TTmcFSQ2ScaleUp.SetTitle("mtb_allhad__ttbar__q2__plus")
	TTmcFSQ2ScaleDown.SetTitle("mtb_allhad__ttbar__q2__minus")

	TTmcFSBup.SetTitle("mtb_allhad__ttbar__btag__plus")
	TTmcFSBdown.SetTitle("mtb_allhad__ttbar__btag__minus")

	TTmcFSTup.SetTitle("mtb_allhad__ttbar__ttag__plus")
	TTmcFSTdown.SetTitle("mtb_allhad__ttbar__ttag__minus")
	#TTmcptup.SetTitle("mtb_allhad__ttbar__pt__plus")
	#TTmcptdown.SetTitle("mtb_allhad__ttbar__pt__minus")

	TTmcFSPtSmearUp.SetTitle("mtb_allhad__ttbar__jer__plus")
	TTmcFSPtSmearDown.SetTitle("mtb_allhad__ttbar__jer__minus")


	TTmcFSPileUp.SetTitle("mtb_allhad__ttbar__pile__plus")
	TTmcFSPileDown.SetTitle("mtb_allhad__ttbar__pile__minus")


	TTmcFSPDFUp.SetTitle("mtb_allhad__ttbar__pdf__plus")
	TTmcFSPDFDown.SetTitle("mtb_allhad__ttbar__pdf__minus")


	TTmcFSTriggerUp.SetTitle("mtb_allhad__ttbar__trig__plus")
	TTmcFSTriggerDown.SetTitle("mtb_allhad__ttbar__trig__minus")

	#STsmcFS.SetTitle("mtb_allhad__sts")
	#STtmcFS.SetTitle("mtb_allhad__stt")
	#STtWmcFS.SetTitle("mtb_allhad__sttw")

        output.cd()

	DataFS.Write("mtb_allhad__DATA")
	DataQCD.Write("mtb_allhad__qcd")
	DataQCDUp.Write("mtb_allhad__qcd__Fit__plus")
	DataQCDDown.Write("mtb_allhad__qcd__Fit__minus")
	DataQCDE1Up.Write("mtb_allhad__qcd__TwoD__plus")
	DataQCDE1Down.Write("mtb_allhad__qcd__TwoD__minus")
	DataQCDE2Up.Write("mtb_allhad__qcd__Alt__plus")
	DataQCDE2Down.Write("mtb_allhad__qcd__Alt__minus")
	DataQCDmodmup.Write("mtb_allhad__qcd__modm__plus")
	DataQCDmodmdown.Write("mtb_allhad__qcd__modm__minus")
	#DataQCDBEH.Write("mtb_allhad__qcd__bkg__plus")
	#DataQCDBEL.Write("mtb_allhad__qcd__bkg__minus")

	TTmcFS.Write("mtb_allhad__ttbar")
	STmcFS.Write("mtb_allhad__st")

	TTmcFSPtScaleUp.Write("mtb_allhad__ttbar__jes__plus")
	TTmcFSPtScaleDown.Write("mtb_allhad__ttbar__jes__minus")

	TTmcFSQ2ScaleUp.Write("mtb_allhad__ttbar__q2__plus")
	TTmcFSQ2ScaleDown.Write("mtb_allhad__ttbar__q2__minus")

	TTmcFSBup.Write("mtb_allhad__ttbar__btag__plus")
	TTmcFSBdown.Write("mtb_allhad__ttbar__btag__minus")

	TTmcFSTup.Write("mtb_allhad__ttbar__ttag__plus")
	TTmcFSTdown.Write("mtb_allhad__ttbar__ttag__minus")

	#TTmcptup.Write("mtb_allhad__ttbar__pt__plus")
	#TTmcptdown.Write("mtb_allhad__ttbar__pt__minus")

	TTmcFSPtSmearUp.Write("mtb_allhad__ttbar__jer__plus")
	TTmcFSPtSmearDown.Write("mtb_allhad__ttbar__jer__minus")

	TTmcFSPileUp.Write("mtb_allhad__ttbar__pile__plus")
	TTmcFSPileDown.Write("mtb_allhad__ttbar__pile__minus")

	TTmcFSPDFUp.Write("mtb_allhad__ttbar__pdf__plus")
	TTmcFSPDFDown.Write("mtb_allhad__ttbar__pdf__minus")


	TTmcFSTriggerUp.Write("mtb_allhad__ttbar__trig__plus")
	TTmcFSTriggerDown.Write("mtb_allhad__ttbar__trig__minus")

	#STsmcFS.Write("mtb_allhad__sts")
	#STtmcFS.Write("mtb_allhad__stt")
	#STtWmcFS.Write("mtb_allhad__sttw")

	for RA in range(0,len(mass)):
		coup = "right"
		SignalB11 = ROOT.TFile("rootfiles/TBanalyzerweightedsignal"+coup+str(mass[RA])+"_Trigger_nominal_none_PSET_"+cuts+".root")

		SignalB11PtScaleUp = ROOT.TFile("rootfiles/TBanalyzerweightedsignal"+coup+str(mass[RA])+"_Trigger_nominal_JES_up_PSET_"+options.cuts+".root")
		SignalB11PtScaleDown = ROOT.TFile("rootfiles/TBanalyzerweightedsignal"+coup+str(mass[RA])+"_Trigger_nominal_JES_down_PSET_"+options.cuts+".root")


		SignalB11PtSmearUp = ROOT.TFile("rootfiles/TBanalyzerweightedsignal"+coup+str(mass[RA])+"_Trigger_nominal_JER_up_PSET_"+options.cuts+".root")
		SignalB11PtSmearDown = ROOT.TFile("rootfiles/TBanalyzerweightedsignal"+coup+str(mass[RA])+"_Trigger_nominal_JER_down_PSET_"+options.cuts+".root")


		SignalB11PileUp = ROOT.TFile("rootfiles/TBanalyzerweightedsignal"+coup+str(mass[RA])+"_Trigger_nominal__pileup_up_PSET_"+options.cuts+".root")
		SignalB11PileDown = ROOT.TFile("rootfiles/TBanalyzerweightedsignal"+coup+str(mass[RA])+"_Trigger_nominal__pileup_down_PSET_"+options.cuts+".root")


		SignalB11PDFUp = ROOT.TFile("rootfiles/TBanalyzerweightedsignal"+coup+str(mass[RA])+"_Trigger_nominal_none_pdf_NNPDF_up_PSET_"+options.cuts+".root")
		SignalB11PDFDown = ROOT.TFile("rootfiles/TBanalyzerweightedsignal"+coup+str(mass[RA])+"_Trigger_nominal_none_pdf_NNPDF_down_PSET_"+options.cuts+".root")



		SignalFS = SignalB11.Get("Mtb")
		SignalBup = SignalB11.Get("MtbBup")
		SignalBdown = SignalB11.Get("MtbBdown")
		SignalTup = SignalB11.Get("MtbTup")
		SignalTdown = SignalB11.Get("MtbTdown")

		SignalTriggerup = SignalB11.Get("Mtbtrigup")
		SignalTriggerdown = SignalB11.Get("Mtbtrigdown")

		SignalPtScaleup = SignalB11PtScaleUp.Get("Mtb")
		SignalPtScaledown = SignalB11PtScaleDown.Get("Mtb")

		SignalPtSmearup = SignalB11PtSmearUp.Get("Mtb")
		SignalPtSmeardown = SignalB11PtSmearDown.Get("Mtb")

		SignalPileup = SignalB11PileUp.Get("Mtb")
		SignalPiledown = SignalB11PileDown.Get("Mtb")

		SignalPDFup = SignalB11PDFUp.Get("Mtb")
		SignalPDFdown = SignalB11PDFDown.Get("Mtb")

		output.cd()

		SignalFS.Rebin(rebin)
		SignalBup.Rebin(rebin)
		SignalBdown.Rebin(rebin)
		SignalTup.Rebin(rebin)
		SignalTdown.Rebin(rebin)
		SignalTriggerup.Rebin(rebin)
		SignalTriggerdown.Rebin(rebin)
		SignalPtScaleup.Rebin(rebin)
		SignalPtScaledown.Rebin(rebin)

		SignalPtSmearup.Rebin(rebin)
		SignalPtSmeardown.Rebin(rebin)


		SignalPileup.Rebin(rebin)
		SignalPiledown.Rebin(rebin)

		SignalPDFup.Rebin(rebin)
		SignalPDFdown.Rebin(rebin)


		SignalFS.SetTitle("mtb_allhad__wp"+str(mass[RA]))
		SignalBup.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__btag__plus")
		SignalBdown.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__btag__minus")
		SignalTup.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__ttag__plus")
		SignalTdown.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__ttag__minus")
		SignalTriggerup.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__trig__plus")
		SignalTriggerdown.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__trig__minus")
		SignalPtScaleup.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__jes__plus")
		SignalPtScaledown.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__jes__minus")
		SignalPtSmearup.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__jer__plus")
		SignalPtSmeardown.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__jer__minus")
		SignalPileup.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__pile__plus")
		SignalPiledown.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__pile__minus")
		SignalPDFup.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__pdf__plus")
		SignalPDFdown.SetTitle("mtb_allhad__wp"+str(mass[RA])+"__pdf__minus")

		SignalFS.SetName("mtb_allhad__wp"+str(mass[RA]))
		SignalBup.SetName("mtb_allhad__wp"+str(mass[RA])+"__btag__plus")
		SignalBdown.SetName("mtb_allhad__wp"+str(mass[RA])+"__btag__minus")
		SignalTup.SetName("mtb_allhad__wp"+str(mass[RA])+"__ttag__plus")
		SignalTdown.SetName("mtb_allhad__wp"+str(mass[RA])+"__ttag__minus")
		SignalTriggerup.SetName("mtb_allhad__wp"+str(mass[RA])+"__trig__plus")
		SignalTriggerdown.SetName("mtb_allhad__wp"+str(mass[RA])+"__trig__minus")
		SignalPtScaleup.SetName("mtb_allhad__wp"+str(mass[RA])+"__jes__plus")
		SignalPtScaledown.SetName("mtb_allhad__wp"+str(mass[RA])+"__jes__minus")
		SignalPtSmearup.SetName("mtb_allhad__wp"+str(mass[RA])+"__jer__plus")
		SignalPtSmeardown.SetName("mtb_allhad__wp"+str(mass[RA])+"__jer__minus")
		SignalPileup.SetName("mtb_allhad__wp"+str(mass[RA])+"__pile__plus")
		SignalPiledown.SetName("mtb_allhad__wp"+str(mass[RA])+"__pile__minus")
		SignalPDFup.SetName("mtb_allhad__wp"+str(mass[RA])+"__pdf__plus")
		SignalPDFdown.SetName("mtb_allhad__wp"+str(mass[RA])+"__pdf__minus")


		SignalFS.Write("mtb_allhad__wp"+str(mass[RA]))
		SignalBup.Write("mtb_allhad__wp"+str(mass[RA])+"__btag__plus")
		SignalBdown.Write("mtb_allhad__wp"+str(mass[RA])+"__btag__minus")
		SignalTup.Write("mtb_allhad__wp"+str(mass[RA])+"__ttag__plus")
		SignalTdown.Write("mtb_allhad__wp"+str(mass[RA])+"__ttag__minus")
		SignalTriggerup.Write("mtb_allhad__wp"+str(mass[RA])+"__trig__plus")
		SignalTriggerdown.Write("mtb_allhad__wp"+str(mass[RA])+"__trig__minus")
		SignalPtScaleup.Write("mtb_allhad__wp"+str(mass[RA])+"__jes__plus")
		SignalPtScaledown.Write("mtb_allhad__wp"+str(mass[RA])+"__jes__minus")
		SignalPtSmearup.Write("mtb_allhad__wp"+str(mass[RA])+"__jer__plus")
		SignalPtSmeardown.Write("mtb_allhad__wp"+str(mass[RA])+"__jer__minus")
		SignalPileup.Write("mtb_allhad__wp"+str(mass[RA])+"__pile__plus")
		SignalPiledown.Write("mtb_allhad__wp"+str(mass[RA])+"__pile__minus")
		SignalPDFup.Write("mtb_allhad__wp"+str(mass[RA])+"__pdf__plus")
		SignalPDFdown.Write("mtb_allhad__wp"+str(mass[RA])+"__pdf__minus")


	
