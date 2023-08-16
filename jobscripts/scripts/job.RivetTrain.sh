#!/bin/bash

AuxParam=("${@:4}")
read -ra AnalysesParam <<< "$AuxParam"

NEvents=${AnalysesParam[0]}
Beam1=${AnalysesParam[1]}
Beam2=${AnalysesParam[2]}
CMSEnergy=${AnalysesParam[3]}
PTHardMin=${AnalysesParam[4]}
PTHardMax=${AnalysesParam[5]}
OutputDir=${AnalysesParam[6]}

outdir="$OutputDir/$SLURM_JOB_ID"
mkdir $outdir

AnalysisDir=${AnalysesParam[9]}
CentralityFile=${AnalysesParam[7]}
CentralityDir=${AnalysesParam[8]}

CentralityPathAndFile="$CentralityDir/$CentralityFile"

if [[ "$Beam1" == "p" ]]; then
  Beam1PDI="2212"
elif [[ "$Beam1" == "Au" ]]; then
  Beam1PDI="1000791970"
elif [[ "$Beam1" == "Cu" ]]; then
  Beam1PDI="1000290630"
elif [[ "$Beam1" == "U" ]]; then
  Beam1PDI="1000922380"
elif [[ "$Beam1" == "He" ]]; then
  Beam1PDI="1000020040"
elif [[ "$Beam1" == "d" ]]; then
  Beam1PDI="1000010020"
elif [[ "$Beam1" == "Pb" ]]; then
  Beam1PDI="1000822080"
fi

if [[ "$Beam2" == "p" ]]; then
  Beam2PDI="2212"
elif [[ "$Beam2" == "Au" ]]; then
  Beam2PDI="1000791970"
elif [[ "$Beam2" == "Cu" ]]; then
  Beam2PDI="1000290630"
elif [[ "$Beam2" == "U" ]]; then
  Beam2PDI="1000922380"
elif [[ "$Beam2" == "He" ]]; then
  Beam2PDI="1000020040"
elif [[ "$Beam2" == "d" ]]; then
  Beam2PDI="1000010020"
elif [[ "$Beam2" == "Pb" ]]; then
  Beam2PDI="1000822080"
fi

#To have completely unique seeds, need to multiply by a number which is greater than the number of jobs.  Protecting in case multiple jobs start at once.
JOBID=`echo $SLURM_JOB_ID`
TASKID=`echo $SLURM_ARRAY_TASK_ID `
NEWTASKNUM=$((TASKID * 10000))
JOB_NUMBER=$((NEWTASKNUM + JOBID))

echo $JOB_NUMBER
mkdir $outdir/$JOB_NUMBER

OUTFILE="RivetFIFO.yoda"

AuxNames=("$2")
AuxFlags=("$3")
AnalysesNames=()
AnalysesFlags=()

for i in ${AuxNames[@]}; do
	AnalysesNames+=("$i")	
done

for i in ${AuxFlags[@]}; do
        AnalysesFlags+=("$i")   
done

source /lustre/isaac/proj/UTK0019/Rivet3.1.8/local/rivetenv.sh

cp /lustre/isaac/proj/UTK0019/sharr100/jobscripts/pythia/pythia8/main93 $outdir/$JOB_NUMBER   
cp /lustre/isaac/proj/UTK0019/sharr100/jobscripts/pythia/pythia8/main93.cmnd $outdir/$JOB_NUMBER 

cp $CentralityPathAndFile $outdir/$JOB_NUMBER
for i in ${!AnalysesNames[@]}; do
	cp $AnalysisDir/${AnalysesNames[$i]}/${AnalysesNames[$i]}.cc $outdir/$JOB_NUMBER
	cp $AnalysisDir/${AnalysesNames[$i]}/${AnalysesNames[$i]}.info $outdir/$JOB_NUMBER
	cp $AnalysisDir/${AnalysesNames[$i]}/${AnalysesNames[$i]}.plot $outdir/$JOB_NUMBER
	cp $AnalysisDir/${AnalysesNames[$i]}/${AnalysesNames[$i]}.yoda $outdir/$JOB_NUMBER
	cp $AnalysisDir/${AnalysesNames[$i]}/Rivet${AnalysesNames[$i]}.so $outdir/$JOB_NUMBER
done

cd $outdir/$JOB_NUMBER

FILES=""

FIFOFILE="$JOB_NUMBER"_"$NEvents"_"$Beam1"_"$Beam2"_"$CMSEnergy".hepmc
FILENAME="$JOB_NUMBER"_"$NEvents"_"$Beam1"_"$Beam2"_"$CMSEnergy"

ANALYSES=""
SEDANALYSES=""
for i in ${!AnalysesNames[@]}; do
	ANALYSES="${AnalysesNames[$i]}${AnalysesFlags[$i]}"
  SEDANALYSES+="${AnalysesNames[$i]}${AnalysesFlags[$i]},"
done
SEDANALYSES=${SEDANALYSES%,}

echo $ANALYSES

# sed -i "s/idA = 2212/idA = $Beam1PDI/g" "main93.cmnd"
# sed -i "s/first beam p = 2212/first beam p = $Beam1PDI/g" "main93.cmnd"
# sed -i "s/idB = 2212/idB = $Beam2PDI/g" "main93.cmnd"
# sed -i "s/second beam p = 2212/second beam p = $Beam2PDI/g" "main93.cmnd"
# sed -i "s/numberOfEvents = 200/numberOfEvents = $NEvents/g" "main93.cmnd"
# sed -i "s/Beams:eCM = 200/Beams:eCM = $CMSEnergy/g" "main93.cmnd"
# sed -i "s/PHENIX_2008_I777211/$SEDANALYSES/g" "main93.cmnd"
# sed -i "s/calibration_PHENIX_AuAu200GeV.yoda/$CentralityFile/g" "main93.cmnd"

#Removes the calibration file for pp collisions
# if [[ "$Beam1" == "p" && "$Beam2" == "p"]]; then
# sed -i '/Main:preload = calibration_PHENIX_AuAu200GeV.yoda/d' main93.cmnd
# fi

#mkfifo /lustre/isaac/proj/UTK0019/sharr100/tmp/$FIFOFILE

./main93 -c main93.cmnd -o $FILENAME -s $JOB_NUMBER  

#rm $FIFOFILE

for i in ${!AnalysesNames[@]}; do
       	rm $outdir/$JOB_NUMBER/${AnalysesNames[$i]}.cc
       	rm $outdir/$JOB_NUMBER/${AnalysesNames[$i]}.info
       	rm $outdir/$JOB_NUMBER/${AnalysesNames[$i]}.plot
       	rm $outdir/$JOB_NUMBER/${AnalysesNames[$i]}.yoda
       	rm $outdir/$JOB_NUMBER/Rivet${AnalysesNames[$i]}.so
		    rm $outdir/$JOB_NUMBER/main93
		    #rm $outdir/$JOB_NUMBER/main93.cmnd
        rm $outdir/$JOB_NUMBER/$CentralityFile
        mv $FILENAME.yoda Rivet.yoda
        rm $outdir/$JOB_NUMBER/$FIFOFILE
done

