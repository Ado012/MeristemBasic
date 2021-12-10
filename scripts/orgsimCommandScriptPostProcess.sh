#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1 # asking for 1 cpus
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=0-02:00:00     # 0 day and 2 hours
#SBATCH --output=my.stdout
#SBATCH --mail-user=ado012@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="orgsim"
#SBATCH -p intel


#MODEL INIT OUTPUT

#./build-OrgSimMerc_IDE-Desktop-Debug/simulator $1 $2 rk6 > $3

mv $1 ./Results/Raw_DataFiles/

Rscript --vanilla ./csvSplitterCommandLine.R "./Results/Raw_DataFiles/results$2" "./Results/Raw_DataFiles/wtresults$2"

python parviewScriptEdit.py $2 "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11"

pvpython weitaoextended7Temp.py



 



