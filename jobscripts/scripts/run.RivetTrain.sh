#!/bin/bash
#SBATCH -J Pythia			       #The name of the job
#SBATCH -A ACF-UTK0019            # The project account to be charged
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1         # cpus per node 
## SBATCH --partition=condo-cnattras
#SBATCH --partition=campus      
#SBATCH --time=0-09:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --output=/lustre/isaac/scratch/sharr100/FIFO/pythia/output/Pythia.o%J
#SBATCH --error=/lustre/isaac/scratch/sharr100/FIFO/pythia/error/Pythia.e%J	 
## SBATCH --qos=condo
#SBATCH --qos=campus

## OPTIONS:
## -A ACF-UTK0019
## --partition=condo-cnattras
## --partition=campus  
## --qos=condo

array=("$2")

for i in ${array[@]}; do
	echo $i
done

srun bash -l -c -v "$SLURM_SUBMIT_DIR/job.RivetTrain.sh \"$1\" \"$2\" \"$3\" \"$4\""