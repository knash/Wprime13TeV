



###################################################################
##								 ##
## Name: Tagrate_Maker_B.py				         ##
## Author: Kevin Nash 						 ##
## Date: 6/5/2012						 ##
## Purpose: This program takes the root files created by  	 ##
##          TBrate.py and creates the average b-tagging rate, 	 ##
##	    then fits the average b-tagging rates		 ##
##          tagrates with a several functions 			 ##
##	    which are stored in the fitdata folder to be used    ##
##	    to weight the pre b tagged sample and create	 ##
##	    QCD background estimates				 ##
##								 ##
###################################################################

import os
import array
import glob
import math
import ROOT
import sys
from array import *
from ROOT import *
# Create the outfiles that will store fit and sigma data for use later
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-s', '--set', metavar='F', type='string', action='store',
                  default	=	'data',
                  dest		=	'set',
                  help		=	'data or QCD')

parser.add_option('-c', '--cuts', metavar='F', type='string', action='store',
                  default	=	'rate',
                  dest		=	'cuts',
                  help		=	'Cuts type (ie default, rate, etc)')

(options, args) = parser.parse_args()

import Wprime_Functions	
from Wprime_Functions import *

rootdir="rootfiles/"

setstr = ""
if options.set=='QCD':
	setstr = 'QCD'


saveout = sys.stdout
Outf1   =   open("fitdata/pol2input"+setstr+"eta1_PSET_"+options.cuts+".txt", "w")
Outf2   =   open("fitdata/pol2input"+setstr+"eta2_PSET_"+options.cuts+".txt", "w")
Outf3   =   open("fitdata/pol2input"+setstr+"eta3_PSET_"+options.cuts+".txt", "w")
Outf4   =   open("fitdata/pol4input"+setstr+"eta1_PSET_"+options.cuts+".txt", "w")
Outf5   =   open("fitdata/pol4input"+setstr+"eta2_PSET_"+options.cuts+".txt", "w")
Outf6   =   open("fitdata/pol4input"+setstr+"eta3_PSET_"+options.cuts+".txt", "w")
Outf7   =   open("fitdata/pol0input"+setstr+"eta1_PSET_"+options.cuts+".txt", "w")
Outf8   =   open("fitdata/pol0input"+setstr+"eta2_PSET_"+options.cuts+".txt", "w")
Outf9   =   open("fitdata/pol0input"+setstr+"eta3_PSET_"+options.cuts+".txt", "w")
Outf10   =   open("fitdata/newfitinput"+setstr+"eta1_PSET_"+options.cuts+".txt", "w")
Outf11   =   open("fitdata/newfitinput"+setstr+"eta2_PSET_"+options.cuts+".txt", "w")
Outf12   =   open("fitdata/newfitinput"+setstr+"eta3_PSET_"+options.cuts+".txt", "w")
Outf13   =   open("fitdata/newfiterrorinput"+setstr+"eta1_PSET_"+options.cuts+".txt", "w")
Outf14   =   open("fitdata/newfiterrorinput"+setstr+"eta2_PSET_"+options.cuts+".txt", "w")
Outf15   =   open("fitdata/newfiterrorinput"+setstr+"eta3_PSET_"+options.cuts+".txt", "w")
Outf16   =   open("fitdata/bpinput"+setstr+"eta1_PSET_"+options.cuts+".txt", "w")
Outf17   =   open("fitdata/bpinput"+setstr+"eta2_PSET_"+options.cuts+".txt", "w")
Outf18   =   open("fitdata/bpinput"+setstr+"eta3_PSET_"+options.cuts+".txt", "w")
Outf19   =   open("fitdata/bperrorinput"+setstr+"eta1_PSET_"+options.cuts+".txt", "w")
Outf20   =   open("fitdata/bperrorinput"+setstr+"eta2_PSET_"+options.cuts+".txt", "w")
Outf21   =   open("fitdata/bperrorinput"+setstr+"eta3_PSET_"+options.cuts+".txt", "w")
Outf22   =   open("fitdata/pol3input"+setstr+"eta1_PSET_"+options.cuts+".txt", "w")
Outf23   =   open("fitdata/pol3input"+setstr+"eta2_PSET_"+options.cuts+".txt", "w")
Outf24   =   open("fitdata/pol3input"+setstr+"eta3_PSET_"+options.cuts+".txt", "w")
Outf25   =   open("fitdata/expoconinput"+setstr+"eta1_PSET_"+options.cuts+".txt", "w")
Outf26   =   open("fitdata/expoconinput"+setstr+"eta2_PSET_"+options.cuts+".txt", "w")
Outf27   =   open("fitdata/expoconinput"+setstr+"eta3_PSET_"+options.cuts+".txt", "w")
Outf28   =   open("fitdata/expoconerrorinput"+setstr+"eta1_PSET_"+options.cuts+".txt", "w")
Outf29   =   open("fitdata/expoconerrorinput"+setstr+"eta2_PSET_"+options.cuts+".txt", "w")
Outf30   =   open("fitdata/expoconerrorinput"+setstr+"eta3_PSET_"+options.cuts+".txt", "w")
Outf31   =   open("fitdata/expolininput"+setstr+"eta1_PSET_"+options.cuts+".txt", "w")
Outf32   =   open("fitdata/expolininput"+setstr+"eta2_PSET_"+options.cuts+".txt", "w")
Outf33   =   open("fitdata/expolininput"+setstr+"eta3_PSET_"+options.cuts+".txt", "w")
Outf34   =   open("fitdata/expolinerrorinput"+setstr+"eta1_PSET_"+options.cuts+".txt", "w")
Outf35   =   open("fitdata/expolinerrorinput"+setstr+"eta2_PSET_"+options.cuts+".txt", "w")
Outf36   =   open("fitdata/expolinerrorinput"+setstr+"eta3_PSET_"+options.cuts+".txt", "w")
sto = sys.stdout  



p0 = 0.0
p1 = 0.0
p2 = 0.0
p3 = 0.0
p4 = 0.0

p00 = 0.0
p01= 0.0
p02= 0.0
p03= 0.0
p04 = 0.0
p11= 0.0
p12= 0.0
p13= 0.0
p14 = 0.0
p22= 0.0
p23= 0.0
p24 = 0.0
p33= 0.0
p34 = 0.0
p44 = 0.0


print "Running on "+options.set
f = TFile(rootdir+"TBratefile"+options.set+"_PSET_"+options.cuts+".root")


f1 = TFile(rootdir+"TBratefilettbar_PSET_"+options.cuts+".root")


SigFiles = [
TFile(rootdir+"TBratefileweightedsignal1300_PSET_"+options.cuts+".root"),
TFile(rootdir+"TBratefileweightedsignal1500_PSET_"+options.cuts+".root"),
TFile(rootdir+"TBratefileweightedsignal1700_PSET_"+options.cuts+".root"),
TFile(rootdir+"TBratefileweightedsignal1900_PSET_"+options.cuts+".root"),
TFile(rootdir+"TBratefileweightedsignal2100_PSET_"+options.cuts+".root"),
TFile(rootdir+"TBratefileweightedsignal2300_PSET_"+options.cuts+".root"),
TFile(rootdir+"TBratefileweightedsignal2700_PSET_"+options.cuts+".root"),
]

output = TFile( "TBrate_Maker_"+setstr+"_PSET_"+options.cuts+".root", "recreate" )
output.cd()

# Get numerators and denominators for each eta region

neta1 = f.Get("pteta1")
deta1 = f.Get("pteta1pretag")

neta2 = f.Get("pteta2")
deta2 = f.Get("pteta2pretag")

neta3 = f.Get("pteta3")
deta3 = f.Get("pteta3pretag")


# Ditto for ttbar

ttneta1 = f1.Get("pteta1")
ttdeta1 = f1.Get("pteta1pretag")

ttneta2 = f1.Get("pteta2")
ttdeta2 = f1.Get("pteta2pretag")

ttneta3 = f1.Get("pteta3")
ttdeta3 = f1.Get("pteta3pretag")



ntot1 = ttneta1.Integral() + neta1.Integral()
ntot2 = ttneta2.Integral() + neta2.Integral()
ntot3 = ttneta3.Integral() + neta3.Integral()

dtot1 = ttdeta1.Integral() + deta1.Integral()
dtot2 = ttdeta2.Integral() + deta2.Integral()
dtot3 = ttdeta3.Integral() + deta3.Integral()

print "pretag QCD & " +strf1(deta1.Integral()) + " ($"+ strf(100*deta1.Integral()/dtot1) + "\%$) & " + strf1(deta2.Integral()) + " ($"+ strf(100*deta2.Integral()/dtot2) + "\%$) & " + strf1(deta3.Integral()) + " ($"+ strf(100*deta3.Integral()/dtot3) + "\%$)"
print "tagged QCD & " +strf1(neta1.Integral()) + " ($"+ strf(100*neta1.Integral()/ntot1) + "\%$) & " + strf1(neta2.Integral()) + " ($"+ strf(100*neta2.Integral()/ntot2) + "\%$) & " + strf1(neta3.Integral()) + " ($"+ strf(100*neta3.Integral()/ntot3) + "\%$)"
print "pretag "+r"$\ttbar$ & " +strf1(ttdeta1.Integral()) + " ($"+ strf(100*ttdeta1.Integral()/dtot1) + "\%$) & " + strf1(ttdeta2.Integral()) + " ($"+ strf(100*ttdeta2.Integral()/dtot2) + "\%$) & " + strf1(ttdeta3.Integral()) + " ($"+ strf(100*ttdeta3.Integral()/dtot3) + "\%$)"
print "tagged "+r"$\ttbar$ & " +strf1(ttneta1.Integral()) + " ($"+ strf(100*ttneta1.Integral()/ntot1) + "\%$) & " + strf1(ttneta2.Integral()) + " ($"+ strf(100*ttneta2.Integral()/ntot2) + "\%$) & " + strf1(ttneta3.Integral()) + " ($"+ strf(100*ttneta3.Integral()/ntot3) + "\%$)"
bins=[]
bins= [370,420,450,500,590,720,1300]
bins2=array('d',bins)

neta1r = neta1.Rebin(len(bins2)-1,"neta1r",bins2)
deta1r = deta1.Rebin(len(bins2)-1,"deta1r",bins2)

neta2r = neta2.Rebin(len(bins2)-1,"neta2r",bins2)
deta2r = deta2.Rebin(len(bins2)-1,"deta2r",bins2)

neta3r = neta3.Rebin(len(bins2)-1,"neta3r",bins2)
deta3r = deta3.Rebin(len(bins2)-1,"deta3r",bins2)

ttneta1r = ttneta1.Rebin(len(bins2)-1,"ttneta1r",bins2)
ttdeta1r = ttdeta1.Rebin(len(bins2)-1,"ttdeta1r",bins2)

ttneta2r = ttneta2.Rebin(len(bins2)-1,"ttneta2r",bins2)
ttdeta2r = ttdeta2.Rebin(len(bins2)-1,"ttdeta2r",bins2)

ttneta3r = ttneta3.Rebin(len(bins2)-1,"ttneta3r",bins2)
ttdeta3r = ttdeta3.Rebin(len(bins2)-1,"ttdeta3r",bins2)

if options.set=='data':
	print 'subtracting ttbar'
	neta1r.Add(ttneta1r,-1)
	deta1r.Add(ttdeta1r,-1)
	neta2r.Add(ttneta2r,-1)
	deta2r.Add(ttdeta2r,-1)
	neta3r.Add(ttneta3r,-1)
	deta3r.Add(ttdeta3r,-1)

outputa = TFile( "B_tagging_sigcont"+setstr+".root", "recreate" )
outputa.cd()
mass = [1300,1500,1700,1900,2100,2300,2700]
for ifile in range(0,len(SigFiles)):
	nseta1 = SigFiles[ifile].Get("pteta1")
	dseta1 = SigFiles[ifile].Get("pteta1pretag")
	nseta2 = SigFiles[ifile].Get("pteta2")
	dseta2 = SigFiles[ifile].Get("pteta2pretag")
	nseta3 = SigFiles[ifile].Get("pteta3")
	dseta3 = SigFiles[ifile].Get("pteta3pretag")
	

	nseta1r = nseta1.Rebin(len(bins2)-1,"nseta1r",bins2)
	dseta1r = dseta1.Rebin(len(bins2)-1,"dseta1r",bins2)

	nseta2r = nseta2.Rebin(len(bins2)-1,"nseta2r",bins2)
	dseta2r = dseta2.Rebin(len(bins2)-1,"dseta2r",bins2)

	nseta3r = nseta3.Rebin(len(bins2)-1,"nseta3r",bins2)
	dseta3r = dseta3.Rebin(len(bins2)-1,"dseta3r",bins2)

	print "pretag signal at "+str(mass[ifile])+"$\GeV$ & " +strf1(dseta1.Integral()) + " ($"+ strf(100*dseta1.Integral()/dtot1) + "\%$) & " + strf1(dseta2.Integral()) + " ($"+ strf(100*dseta2.Integral()/dtot2) + "\%$) & " + strf1(dseta3.Integral()) + " ($"+ strf(100*dseta3.Integral()/dtot3) + "\%$)"
	print "tagged signal at "+str(mass[ifile])+"$\GeV$ & " +strf1(nseta1.Integral()) + " ($"+ strf(100*nseta1.Integral()/ntot1) + "\%$) & " + strf1(nseta2.Integral()) + " ($"+ strf(100*nseta2.Integral()/ntot2) + "\%$) & " + strf1(nseta3.Integral()) + " ($"+ strf(100*nseta3.Integral()/ntot3) + "\%$)"

	nseta1r.Add(neta1r)
	dseta1r.Add(deta1r)

	nseta2r.Add(neta2r)
	dseta2r.Add(deta2r)

	nseta3r.Add(neta3r)
	dseta3r.Add(deta3r)

	tagrateseta1 = nseta1r.Clone("tagrateseta1"+str(mass[ifile]))
	tagrateseta1.Divide(tagrateseta1,dseta1r,1.0,1.0,"B")

	tagrateseta2 = nseta2r.Clone("tagrateseta2"+str(mass[ifile]))
	tagrateseta2.Divide(tagrateseta2,dseta2r,1.0,1.0,"B")

	tagrateseta3 = nseta3r.Clone("tagrateseta3"+str(mass[ifile]))
	tagrateseta3.Divide(tagrateseta3,dseta3r,1.0,1.0,"B")


	tagrateseta1.Write()

	tagrateseta2.Write()

	tagrateseta3.Write()








output.cd()

# Create subtracted tagrates by division

tagrateeta1 = neta1r.Clone("tagrateeta1")
#tagrateeta1.Divide(deta1r)
tagrateeta1.Divide(tagrateeta1,deta1r,1.0,1.0,"B")
tagrateeta2 = neta2r.Clone("tagrateeta2")
#tagrateeta2.Divide(deta2r)
tagrateeta2.Divide(tagrateeta2,deta2r,1.0,1.0,"B")
tagrateeta3 = neta3r.Clone("tagrateeta3")
#tagrateeta3.Divide(deta3r)
tagrateeta3.Divide(tagrateeta3,deta3r,1.0,1.0,"B")


tagrateeta1.Write()
tagrateeta2.Write()
tagrateeta3.Write()

# The rest here writes the fit, and uses the covariance matrix to propagate errors for the pol4 function

sys.stdout = saveout
print "------------------------------------"
print "POL2"
print "------------------------------------"

sys.stdout = Outf1
tagrateeta1.Fit("pol2","F")
fitter = TVirtualFitter.GetFitter()

for i in range(0,3):
	print(fitter.GetParameter(i))

p0 = fitter.GetCovarianceMatrixElement(0,0)
p1 = 2*fitter.GetCovarianceMatrixElement(0,1)
p2 = 2*fitter.GetCovarianceMatrixElement(0,2) + fitter.GetCovarianceMatrixElement(1,1)
p3 = 2*fitter.GetCovarianceMatrixElement(2,1)
p4 = fitter.GetCovarianceMatrixElement(2,2)

sys.stdout = Outf4
print str(p0)
print str(p1)
print str(p2)
print str(p3)
print str(p4)

sys.stdout = Outf2
tagrateeta2.Fit("pol2","F")
fitter = TVirtualFitter.GetFitter()
for i in range(0,3):
	print(fitter.GetParameter(i))

p0 = fitter.GetCovarianceMatrixElement(0,0)
p1 = 2*fitter.GetCovarianceMatrixElement(0,1)
p2 = 2*fitter.GetCovarianceMatrixElement(0,2) + fitter.GetCovarianceMatrixElement(1,1)
p3 = 2*fitter.GetCovarianceMatrixElement(2,1)
p4 = fitter.GetCovarianceMatrixElement(2,2)

sys.stdout = Outf5
print str(p0)
print str(p1)
print str(p2)
print str(p3)
print str(p4)

sys.stdout = Outf3
tagrateeta3.Fit("pol2","F")
fitter = TVirtualFitter.GetFitter()
for i in range(0,3):
	print(fitter.GetParameter(i))

p0 = fitter.GetCovarianceMatrixElement(0,0)
p1 = 2*fitter.GetCovarianceMatrixElement(0,1)
p2 = 2*fitter.GetCovarianceMatrixElement(0,2) + fitter.GetCovarianceMatrixElement(1,1)
p3 = 2*fitter.GetCovarianceMatrixElement(2,1)
p4 = fitter.GetCovarianceMatrixElement(2,2)

sys.stdout = Outf6
print str(p0)
print str(p1)
print str(p2)
print str(p3)
print str(p4)
sys.stdout = saveout
print "------------------------------------"
print "POL0"
print "------------------------------------"
sys.stdout = Outf7
tagrateeta1.Fit("pol0","F")
fitter = TVirtualFitter.GetFitter()
print(fitter.GetParameter(0))

sys.stdout = Outf8
tagrateeta2.Fit("pol0","F")
fitter = TVirtualFitter.GetFitter()
print(fitter.GetParameter(0))

sys.stdout = Outf9
tagrateeta3.Fit("pol0","F")
fitter = TVirtualFitter.GetFitter()
print(fitter.GetParameter(0))

sys.stdout = saveout
print "------------------------------------"
print "BIFPOLY"
print "------------------------------------"

BIFP=500.0
BP =TF1("BP",BifPoly,370,1400,5)
BP.FixParameter(4,BIFP)
#fixbin = tagrateeta1.FindBin(BIFP)
#fix = tagrateeta1.GetBinContent(fixbin)
#BP.FixParameter(0,fix)
c4 = TCanvas('c4', 'Tagrate1', 1300, 600)
c4.cd()
tagrateeta1.Fit("BP","F")
#print "CHISQUARE"
#print "chsq"
#print BP.GetChisquare()
#print "ndof"
#print BP.GetNDF()
#print "ratio"
#print BP.GetChisquare()/BP.GetNDF()
#print ""
sys.stdout = Outf16
fitter = TVirtualFitter.GetFitter()
print(fitter.GetParameter(0))
print(fitter.GetParameter(1))
print(fitter.GetParameter(2))
print(fitter.GetParameter(3))
print(BIFP)
print BP.GetChisquare()
print BP.GetNDF()
sys.stdout = Outf19
p00 = fitter.GetCovarianceMatrixElement(0,0)
p01 = 2*fitter.GetCovarianceMatrixElement(0,1)
p02 = 2*fitter.GetCovarianceMatrixElement(0,2)
p03 = 2*fitter.GetCovarianceMatrixElement(0,3)
p11 = fitter.GetCovarianceMatrixElement(1,1)
p12 = 2*fitter.GetCovarianceMatrixElement(1,2)
p13 = 2*fitter.GetCovarianceMatrixElement(1,3)
p22 = fitter.GetCovarianceMatrixElement(2,2)
p23 = 2*fitter.GetCovarianceMatrixElement(2,3)
p33 = fitter.GetCovarianceMatrixElement(3,3)
#The funky order here is important 
print str(p00)
print str(p11)
print str(p22)
print str(p01)
print str(p02)
print str(p12)
print str(p33)
print str(p03)
print str(p13)
print(BIFP)



tagrateeta1.Draw()
c4.Print("BPTAGETA1FIT"+setstr+".root","root")
#fixbin = tagrateeta2.FindBin(BIFP)
#fix = tagrateeta2.GetBinContent(fixbin)
#BP.FixParameter(0,fix)
c3 = TCanvas('c3', 'Tagrate2', 1300, 600)
c3.cd()

BIFP=500.0
BP =TF1("BP",BifPoly,370,1400,5)
BP.FixParameter(4,BIFP)
tagrateeta2.Fit("BP","F")
sys.stdout = saveout
#print "CHISQUARE"
#print "chsq"
#print BP.GetChisquare()
#print "ndof"
#print BP.GetNDF()
#print "ratio"
#print BP.GetChisquare()/BP.GetNDF()
#print ""
sys.stdout = Outf17
fitter = TVirtualFitter.GetFitter()
print(fitter.GetParameter(0))
print(fitter.GetParameter(1))
print(fitter.GetParameter(2))
print(fitter.GetParameter(3))
print(BIFP)
print BP.GetChisquare()
print BP.GetNDF()
sys.stdout = Outf20
p00 = fitter.GetCovarianceMatrixElement(0,0)
p01 = 2*fitter.GetCovarianceMatrixElement(0,1)
p02 = 2*fitter.GetCovarianceMatrixElement(0,2)
p03 = 2*fitter.GetCovarianceMatrixElement(0,3)
p11 = fitter.GetCovarianceMatrixElement(1,1)
p12 = 2*fitter.GetCovarianceMatrixElement(1,2)
p13 = 2*fitter.GetCovarianceMatrixElement(1,3)
p22 = fitter.GetCovarianceMatrixElement(2,2)
p23 = 2*fitter.GetCovarianceMatrixElement(2,3)
p33 = fitter.GetCovarianceMatrixElement(3,3)
print str(p00)
print str(p11)
print str(p22)
print str(p01)
print str(p02)
print str(p12)
print str(p33)
print str(p03)
print str(p13)
print(BIFP)
tagrateeta2.Draw()
c3.Print("BPTAGETA2FIT"+setstr+".root","root")
c2 = TCanvas('c2', 'Tagrate3', 1300, 600)
c2.cd()
#fixbin = tagrateeta3.FindBin(BIFP)
#fix = tagrateeta3.GetBinContent(fixbin)
#BP.FixParameter(0,fix)

BIFP=550.0
BP =TF1("BP",BifPoly,370,1400,5)
BP.FixParameter(4,BIFP)
tagrateeta3.Fit("BP","F")
sys.stdout = saveout
#print "CHISQUARE"
#print "chsq"
#print BP.GetChisquare()
#print "ndof"
#print BP.GetNDF()
#print "ratio"
#print BP.GetChisquare()/BP.GetNDF()
#print ""
sys.stdout = Outf18
fitter = TVirtualFitter.GetFitter()
print(fitter.GetParameter(0))
print(fitter.GetParameter(1))
print(fitter.GetParameter(2))
print(fitter.GetParameter(3))
print(BIFP)
print BP.GetChisquare()
print BP.GetNDF()
sys.stdout = Outf21
p00 = fitter.GetCovarianceMatrixElement(0,0)
p01 = 2*fitter.GetCovarianceMatrixElement(0,1)
p02 = 2*fitter.GetCovarianceMatrixElement(0,2)
p03 = 2*fitter.GetCovarianceMatrixElement(0,3)
p11 = fitter.GetCovarianceMatrixElement(1,1)
p12 = 2*fitter.GetCovarianceMatrixElement(1,2)
p13 = 2*fitter.GetCovarianceMatrixElement(1,3)
p22 = fitter.GetCovarianceMatrixElement(2,2)
p23 = 2*fitter.GetCovarianceMatrixElement(2,3)
p33 = fitter.GetCovarianceMatrixElement(3,3)
print str(p00)
print str(p11)
print str(p22)
print str(p01)
print str(p02)
print str(p12)
print str(p33)
print str(p03)
print str(p13)
print(BIFP)
tagrateeta3.Draw()
c2.Print("BPTAGETA3FIT"+setstr+".root","root")
sys.stdout = saveout
print "------------------------------------"
print "POL3"
print "------------------------------------"

#c1.cd()
sys.stdout = Outf22
tagrateeta1.Fit("pol3","F")
fitter = TVirtualFitter.GetFitter()
for i in range(0,4):
	print(fitter.GetParameter(i))
#tagrateeta1.Draw()

#c2.cd()
sys.stdout = Outf23
tagrateeta2.Fit("pol3","F")
fitter = TVirtualFitter.GetFitter()
for i in range(0,4):
	print(fitter.GetParameter(i))
#tagrateeta2.Draw()

#c3.cd()
sys.stdout = Outf24
tagrateeta3.Fit("pol3","F")
fitter = TVirtualFitter.GetFitter()
for i in range(0,4):
	print(fitter.GetParameter(i))
#tagrateeta3.Draw()


expofitlin =TF1("expofitlin","expo(0) + pol1(2)")
sys.stdout = saveout
print "------------------------------------"
print "EXPO + LINEAR"
print "------------------------------------"

#expofitlin.SetParameter(0, -1.0e+00)
#expofitlin.SetParameter(1, -7.4e-03)
#expofitlin.SetParameter(2, 4.1e-02)
#expofitlin.SetParameter(3, 4.0e-06)

#c1.cd()
sys.stdout = Outf31
tagrateeta1.Fit("expofitlin","F")
fitter = TVirtualFitter.GetFitter()
for i in range(0,4):
	print(fitter.GetParameter(i))
sys.stdout = Outf34
p00 = fitter.GetCovarianceMatrixElement(0,0)
p01 = 2*fitter.GetCovarianceMatrixElement(0,1)
p02 = 2*fitter.GetCovarianceMatrixElement(0,2)
p03 = 2*fitter.GetCovarianceMatrixElement(0,3)
p04 = 2*fitter.GetCovarianceMatrixElement(0,4)
p11 = fitter.GetCovarianceMatrixElement(1,1)
p12 = 2*fitter.GetCovarianceMatrixElement(1,2)
p13 = 2*fitter.GetCovarianceMatrixElement(1,3)
p14 = 2*fitter.GetCovarianceMatrixElement(1,4)
p22 = fitter.GetCovarianceMatrixElement(2,2)
p23 = 2*fitter.GetCovarianceMatrixElement(2,3)
p24 = 2*fitter.GetCovarianceMatrixElement(2,4)
p33 = fitter.GetCovarianceMatrixElement(3,3)
p34 = 2*fitter.GetCovarianceMatrixElement(3,4)
p44 = fitter.GetCovarianceMatrixElement(4,4)
print str(p00)
print str(p01)
print str(p02)
print str(p03)
print str(p04)
print str(p11)
print str(p12)
print str(p13)
print str(p14)
print str(p22)
print str(p23)
print str(p24)
print str(p33)
print str(p34)
print str(p44)
#tagrateeta1.Draw()

#c2.cd()
sys.stdout = Outf32
tagrateeta2.Fit("expofitlin","F")
fitter = TVirtualFitter.GetFitter()
for i in range(0,4):
	print(fitter.GetParameter(i))
sys.stdout = Outf35
p00 = fitter.GetCovarianceMatrixElement(0,0)
p01 = 2*fitter.GetCovarianceMatrixElement(0,1)
p02 = 2*fitter.GetCovarianceMatrixElement(0,2)
p03 = 2*fitter.GetCovarianceMatrixElement(0,3)
p04 = 2*fitter.GetCovarianceMatrixElement(0,4)
p11 = fitter.GetCovarianceMatrixElement(1,1)
p12 = 2*fitter.GetCovarianceMatrixElement(1,2)
p13 = 2*fitter.GetCovarianceMatrixElement(1,3)
p14 = 2*fitter.GetCovarianceMatrixElement(1,4)
p22 = fitter.GetCovarianceMatrixElement(2,2)
p23 = 2*fitter.GetCovarianceMatrixElement(2,3)
p24 = 2*fitter.GetCovarianceMatrixElement(2,4)
p33 = fitter.GetCovarianceMatrixElement(3,3)
p34 = 2*fitter.GetCovarianceMatrixElement(3,4)
p44 = fitter.GetCovarianceMatrixElement(4,4)
print str(p00)
print str(p01)
print str(p02)
print str(p03)
print str(p04)
print str(p11)
print str(p12)
print str(p13)
print str(p14)
print str(p22)
print str(p23)
print str(p24)
print str(p33)
print str(p34)
print str(p44)
#tagrateeta2.Draw()

#c1.cd()
sys.stdout = Outf33
tagrateeta3.Fit("expofitlin","F")
fitter = TVirtualFitter.GetFitter()
for i in range(0,4):
	print(fitter.GetParameter(i))
sys.stdout = Outf36
p00 = fitter.GetCovarianceMatrixElement(0,0)
p01 = 2*fitter.GetCovarianceMatrixElement(0,1)
p02 = 2*fitter.GetCovarianceMatrixElement(0,2)
p03 = 2*fitter.GetCovarianceMatrixElement(0,3)
p04 = 2*fitter.GetCovarianceMatrixElement(0,4)
p11 = fitter.GetCovarianceMatrixElement(1,1)
p12 = 2*fitter.GetCovarianceMatrixElement(1,2)
p13 = 2*fitter.GetCovarianceMatrixElement(1,3)
p14 = 2*fitter.GetCovarianceMatrixElement(1,4)
p22 = fitter.GetCovarianceMatrixElement(2,2)
p23 = 2*fitter.GetCovarianceMatrixElement(2,3)
p24 = 2*fitter.GetCovarianceMatrixElement(2,4)
p33 = fitter.GetCovarianceMatrixElement(3,3)
p34 = 2*fitter.GetCovarianceMatrixElement(3,4)
p44 = fitter.GetCovarianceMatrixElement(4,4)
print str(p00)
print str(p01)
print str(p02)
print str(p03)
print str(p04)
print str(p11)
print str(p12)
print str(p13)
print str(p14)
print str(p22)
print str(p23)
print str(p24)
print str(p33)
print str(p34)
print str(p44)
#tagrateeta3.Draw()
sys.stdout = saveout
print "------------------------------------"
print "EXPO + CONSTANT"
print "------------------------------------"
expofit =TF1("expofit","expo(0) + pol0(2)")

#expofit.SetParameter(0, -1.03e+00);
#expofit.SetParameter(1, 0.05);
#expofit.SetParameter(2, 0.2);

#c1.cd()
sys.stdout = Outf25
tagrateeta1.Fit("expofit","F")
fitter = TVirtualFitter.GetFitter()
for i in range(0,3):
	print(fitter.GetParameter(i))
sys.stdout = Outf28
p00 = fitter.GetCovarianceMatrixElement(0,0)
p01 = 2*fitter.GetCovarianceMatrixElement(0,1)
p02 = 2*fitter.GetCovarianceMatrixElement(0,2)
p03 = 2*fitter.GetCovarianceMatrixElement(0,3)
p11 = fitter.GetCovarianceMatrixElement(1,1)
p12 = 2*fitter.GetCovarianceMatrixElement(1,2)
p13 = 2*fitter.GetCovarianceMatrixElement(1,3)
p22 = fitter.GetCovarianceMatrixElement(2,2)
p23 = 2*fitter.GetCovarianceMatrixElement(2,3)
p33 = fitter.GetCovarianceMatrixElement(3,3)
print str(p00)
print str(p01)
print str(p02)
print str(p03)
print str(p11)
print str(p12)
print str(p13)
print str(p22)
print str(p23)
print str(p33)
#tagrateeta1.Draw()

c2.cd()
sys.stdout = Outf26
tagrateeta2.Fit("expofit","F")
fitter = TVirtualFitter.GetFitter()
for i in range(0,3):
	print(fitter.GetParameter(i))
sys.stdout = Outf29
p00 = fitter.GetCovarianceMatrixElement(0,0)
p01 = 2*fitter.GetCovarianceMatrixElement(0,1)
p02 = 2*fitter.GetCovarianceMatrixElement(0,2)
p03 = 2*fitter.GetCovarianceMatrixElement(0,3)
p11 = fitter.GetCovarianceMatrixElement(1,1)
p12 = 2*fitter.GetCovarianceMatrixElement(1,2)
p13 = 2*fitter.GetCovarianceMatrixElement(1,3)
p22 = fitter.GetCovarianceMatrixElement(2,2)
p23 = 2*fitter.GetCovarianceMatrixElement(2,3)
p33 = fitter.GetCovarianceMatrixElement(3,3)
print str(p00)
print str(p01)
print str(p02)
print str(p03)
print str(p11)
print str(p12)
print str(p13)
print str(p22)
print str(p23)
print str(p33)
#tagrateeta2.Draw()

#c1.cd()
sys.stdout = Outf27
tagrateeta3.Fit("expofit","F")
fitter = TVirtualFitter.GetFitter()
for i in range(0,3):
	print(fitter.GetParameter(i))
sys.stdout = Outf30
p00 = fitter.GetCovarianceMatrixElement(0,0)
p01 = 2*fitter.GetCovarianceMatrixElement(0,1)
p02 = 2*fitter.GetCovarianceMatrixElement(0,2)
p03 = 2*fitter.GetCovarianceMatrixElement(0,3)
p11 = fitter.GetCovarianceMatrixElement(1,1)
p12 = 2*fitter.GetCovarianceMatrixElement(1,2)
p13 = 2*fitter.GetCovarianceMatrixElement(1,3)
p22 = fitter.GetCovarianceMatrixElement(2,2)
p23 = 2*fitter.GetCovarianceMatrixElement(2,3)
p33 = fitter.GetCovarianceMatrixElement(3,3)
print str(p00)
print str(p01)
print str(p02)
print str(p03)
print str(p11)
print str(p12)
print str(p13)
print str(p22)
print str(p23)
print str(p33)
#tagrateeta3.Draw()

sys.stdout = saveout
print "------------------------------------"
print "FIT"
print "------------------------------------"
FIT =TF1("FIT","[0]*([1]+x)/([2]+x)+[3]*x")
sys.stdout = Outf10
tagrateeta1.Fit("FIT","F")
fitter = TVirtualFitter.GetFitter()
print(fitter.GetParameter(0))
print(fitter.GetParameter(1))
print(fitter.GetParameter(2))
print(fitter.GetParameter(3))

p00 = fitter.GetCovarianceMatrixElement(0,0)
p01 = 2*fitter.GetCovarianceMatrixElement(0,1)
p02 = 2*fitter.GetCovarianceMatrixElement(0,2)
p03 = 2*fitter.GetCovarianceMatrixElement(0,3)
p11 = fitter.GetCovarianceMatrixElement(1,1)
p12 = 2*fitter.GetCovarianceMatrixElement(1,2)
p13 = 2*fitter.GetCovarianceMatrixElement(1,3)
p22 = fitter.GetCovarianceMatrixElement(2,2)
p23 = 2*fitter.GetCovarianceMatrixElement(2,3)
p33 = fitter.GetCovarianceMatrixElement(3,3)
sys.stdout = Outf13
#print_options.set_float_precision(2)
print str(p00)
print str(p01)
print str(p02)
print str(p03)
print str(p11)
print str(p12)
print str(p13)
print str(p22)
print str(p23)
print str(p33)

sys.stdout = Outf11
tagrateeta2.Fit("FIT","F")
fitter = TVirtualFitter.GetFitter()
print(fitter.GetParameter(0))
print(fitter.GetParameter(1))
print(fitter.GetParameter(2))
print(fitter.GetParameter(3))

sys.stdout = Outf14
p00 = fitter.GetCovarianceMatrixElement(0,0)
p01 = 2*fitter.GetCovarianceMatrixElement(0,1)
p02 = 2*fitter.GetCovarianceMatrixElement(0,2)
p03 = 2*fitter.GetCovarianceMatrixElement(0,3)
p11 = fitter.GetCovarianceMatrixElement(1,1)
p12 = 2*fitter.GetCovarianceMatrixElement(1,2)
p13 = 2*fitter.GetCovarianceMatrixElement(1,3)
p22 = fitter.GetCovarianceMatrixElement(2,2)
p23 = 2*fitter.GetCovarianceMatrixElement(2,3)
p33 = fitter.GetCovarianceMatrixElement(3,3)
print str(p00)
print str(p01)
print str(p02)
print str(p03)
print str(p11)
print str(p12)
print str(p13)
print str(p22)
print str(p23)
print str(p33)

sys.stdout = Outf12
tagrateeta3.Fit("FIT","F")
fitter = TVirtualFitter.GetFitter()
print(fitter.GetParameter(0))
print(fitter.GetParameter(1))
print(fitter.GetParameter(2))
print(fitter.GetParameter(3))

sys.stdout = Outf15
p00 = fitter.GetCovarianceMatrixElement(0,0)
p01 = 2*fitter.GetCovarianceMatrixElement(0,1)
p02 = 2*fitter.GetCovarianceMatrixElement(0,2)
p03 = 2*fitter.GetCovarianceMatrixElement(0,3)
p11 = fitter.GetCovarianceMatrixElement(1,1)
p12 = 2*fitter.GetCovarianceMatrixElement(1,2)
p13 = 2*fitter.GetCovarianceMatrixElement(1,3)
p22 = fitter.GetCovarianceMatrixElement(2,2)
p23 = 2*fitter.GetCovarianceMatrixElement(2,3)
p33 = fitter.GetCovarianceMatrixElement(3,3)
print str(p00)
print str(p01)
print str(p02)
print str(p03)
print str(p11)
print str(p12)
print str(p13)
print str(p22)
print str(p23)
print str(p33)


