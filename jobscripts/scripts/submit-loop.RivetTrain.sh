#!/bin/bash
#AnalysesNames=("PHENIX_2008_I777211 PHENIX_2009_I816486")
AnalysesNames=("PHENIX_2011_I886590")
#Flag example ":cent=GEN" "->" for pythia8 cmnd file (beam example: AUAU200 or PP200)
#AnalysesFlags=(":cent->GEN:beam->AUAU200 :cent->GEN:beam->AUAU200")
AnalysesFlags=(":beam->PP62.4")
OutputDir="/lustre/isaac/proj/UTK0019/sharr100/pythia_output"
CentralityFile="calibration_PHENIX_PP200GeV.yoda"
CentralityDir="/lustre/isaac/proj/UTK0019/sharr100/rivet/RIVETAnalyses/Centralities/Calibration"
AnalysisDir="/lustre/isaac/proj/UTK0019/sharr100/rivet/RIVETAnalyses"
#Parameters: Number of events, beam1, beam2, cms energy, min pT-hard, max pT-hard
AnalysesParameters=("5000 p p 62.4 0 -1 $OutputDir $CentralityFile $CentralityDir $AnalysisDir")
for (( i=1; i<=1; i++ ))
do
   source submit.RivetTrain.sh $i "${AnalysesNames[*]}" "${AnalysesFlags[*]}" "${AnalysesParameters[*]}"
done
