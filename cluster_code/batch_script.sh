#PBS -l nodes=1:ppn=1,walltime=23:00:00
#PBS -N aggregates_algorithms_testing_Worm1_001
#PBS -M plutzer@wuslt.edu
#PBS -d /home/plutzer/
#PBS -e /home/plutzer/job_log_files/
#PBS -o /home/plutzer/job_log_files/
cd /home/plutzer/code

export ANACONDA_PATH="/act/Anaconda3-2.3.0/bin"
export PATH=$ANACONDA_PATH:$PATH

source activate zplab_cluster

python aggregates_crucible.py
