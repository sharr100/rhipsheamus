#!/bin/bash
AnalysisDir="/lustre/isaac/proj/UTK0019/sharr100/rivet/RIVETAnalyses"
AnalysesNames=("PHENIX_2007_I731133")
#AnalysesNames="PHENIX_2003_I619987"
#OUTDIR="/lustre/isaac/proj/UTK0019/sharr100/pythia_output"
OUTDIR="/lustre/isaac/proj/UTK0019/sharr100/pythia_output"

for i in ${!AnalysesNames[@]}; do

    cd "$AnalysisDir/${AnalysesNames[$i]}/"
    #rivet-buildplugin "Rivet${AnalysesNames[$i]}.so" "${AnalysesNames[$i]}.cc"

    FILES=""
    for d in $OUTDIR/*/* ; do
        if test -f "$d/Rivet.yoda"; then
    	    FILES="$FILES $d/Rivet.yoda"
        fi
    done
    
    echo $FILES
    rivet-merge --pwd -o Angantyr.yoda -O beam $FILES
    rivet-mkhtml --pwd Angantyr.yoda --errs
    #rivet-mkhtml --pwd Rivet_final.yoda --match "d\d\d-x\d\d-y\d\d$"
done
