	#PBS -l nodes=1:ppn=1,walltime=47:00:00
	#PBS -N aggregates_algorithms_testing_Worm1_001
	#PBS -M plutzer@wuslt.edu
	#PBS -d /scratch/plutzer/
	#PBS -e /home/plutzer/job_log_files/
	#PBS -o /home/plutzer/job_log_files/
	# #PBS -q old

	export ANACONDA_PATH="/act/Anaconda3-2.3.0/bin"
	export PATH=$ANACONDA_PATH:$PATH

	source activate zplab_cluster
	cd $PBS_O_WORKDIR
	python /home/plutzer/code/aggregates_crucible.py
