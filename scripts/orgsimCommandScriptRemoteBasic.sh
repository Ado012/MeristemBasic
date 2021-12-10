#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1 # asking for 1 cpus
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-03:00:00     # 0 day and 3 hours
#SBATCH --output=my.stdout
#SBATCH --mail-user=ado012@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="orgsim"
#SBATCH -p intel


#MODEL INIT OUTPUT

./build-OrgSimBasic-Desktop-Debug/simulator $1 $2 rk6 > $3

mv $3 ./Results/Raw_DataFiles/

Rscript --vanilla ./csvSplitterCommandLine.R "./Results/Raw_DataFiles/results$4" "./Results/Raw_DataFiles/wtresults$4"

python parviewScriptEditRemote.py "wtresults$4"

./ParaView/bin/pvpython weitaoextendedRemote.py




 



