
cd /data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG

## activate conda env
# conda activate atac
sinteractive \
 --account=punim0586  \
 --ntasks=1 \
 --cpus-per-task=1 \
 --mem=80000 \
 --time=8:00:00 \
 --partition=mig
 
# list of modules
source /usr/local/module/spartan_new.sh
module load web_proxy
module load gcc/8.3.0 openmpi/3.1.4
module load r/4.0.0  
