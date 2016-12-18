#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

ch=(ele mu)
samplePath=/data7/htong/skim_NCUGlobalTuples

for ((i=0; i<${#ch[@]}; i++)); do

    if [ `echo ${ch[$i]} | grep -c "mu"` -gt 0 ]; then

	echo "Processing muon data set..."

	root -q -b -l alphaTree.C\(\"$samplePath/skim_mu_SingleMuon_Run2015D-05Oct2015-v1.root\"\,\"SingleMuon-Run2015D-v1\"\,\"mu\"\)
	root -q -b -l alphaTree.C\(\"$samplePath/skim_mu_SingleMuon_Run2015D-PromptReco-v4.root\"\,\"SingleMuon-Run2015D-v4\"\,\"mu\"\)

    elif [ `echo ${ch[$i]} | grep -c "ele"` -gt 0 ]; then

	echo "Processing electron data set..."

        root -q -b -l alphaTree.C\(\"$samplePath/skim_ele_SingleElectron_Run2015D-05Oct2015-v1.root\"\,\"SingleElectron-Run2015D-v1\"\,\"ele\"\)
	root -q -b -l alphaTree.C\(\"$samplePath/skim_ele_SingleElectron_Run2015D-PromptReco-v4.root\"\,\"SingleElectron-Run2015D-v4\"\,\"ele\"\)

    fi

    mv *root data

    echo "Processing Z+jets background..."

    root -q -b -l alphaTree.C\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-100to200_13TeV\"\,\"${ch[$i]}\"\)
    root -q -b -l alphaTree.C\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-200to400_13TeV\"\,\"${ch[$i]}\"\)
    root -q -b -l alphaTree.C\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-400to600_13TeV\"\,\"${ch[$i]}\"\)
    root -q -b -l alphaTree.C\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-600toInf_13TeV\"\,\"${ch[$i]}\"\)

    mv *root Zjets

    echo "Processing VV background..."

    root -q -b -l alphaTree.C\(\"$samplePath/skim_${ch[$i]}_crab_WW_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WW_TuneCUETP8M1_13TeV\"\,\"${ch[$i]}\"\)
    root -q -b -l alphaTree.C\(\"$samplePath/skim_${ch[$i]}_crab_WZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WZ_TuneCUETP8M1_13TeV\"\,\"${ch[$i]}\"\)
    root -q -b -l alphaTree.C\(\"$samplePath/skim_${ch[$i]}_crab_ZZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"ZZ_TuneCUETP8M1_13TeV\"\,\"${ch[$i]}\"\)

    echo "Processing ZH background..."

    root -q -b -l alphaTree.C\(\"$samplePath/skim_${ch[$i]}_crab_ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8.root\"\,\"ZH_HToBB_ZToLL_M125_13TeV\"\,\"${ch[$i]}\"\)

    echo "Processing ttbar background..."

    root -q -b -l alphaTree.C\(\"$samplePath/skim_${ch[$i]}_crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root\"\,\"TT_TuneCUETP8M1_13TeV\"\,\"${ch[$i]}\"\)

    mv *root minor

    rm -f inputdir.txt
    rm -f *.pcm *.d *.so

    echo "Done."

done

echo "All the jobs are finished."

exit