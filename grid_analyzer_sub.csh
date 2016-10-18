#! /bin/sh
tar czvf tarball.tgz fitdata Tagrate2D*.root bkgwQCD*.root ModMassFile_PSET_rate*.root rootlogon.C TBanalyzer.py TBkinematics.py Wprime_Functions.py TBrate.py  Triggerweight_databtagsmass.root Triggerweight_signalright1200btags.root PileUp_Ratio_ttbar.root PileUp_Ratio_signal*.root 
./development/runManySections.py --createCommandFile --cmssw --addLog --setTarball=tarball.tgz \ana.listOfJobs commands.cmd
./runManySections.py --submitCondor commands.cmd
condor_q knash
