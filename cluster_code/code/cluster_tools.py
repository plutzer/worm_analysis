import json
import os
import pathlib
import subprocess
import string
import sys
sys.path.append('/home/plutzer/code/wormPhysiology/measurePhysiology/')
import wormPhysiology.measurePhysiology.organizeData as organizeData
import IsaacMeasure_Health as hm
#import wormPhysiology.measurePhysiology.IsaacMeasure_Health as hm

# Constants
file_dir = '/scratch/plutzer/'   # Where files will be output for running jobs
work_dir = '/scratch/plutzer/20180628_Q40_4/work_dir'  # Where mask data will be saved
human_dir = '/home/plutzer/utilities/' # Where data for texture classifiers are stored

def configure_pipeline(*dir_args):
    # dir_args - set of Path variables
    animals_toprocess = []
    
    print("Enumerating experiment directories to process")
    for data_dir in dir_args:
        checker = organizeData.HumanCheckpoints(data_dir)
        checker.clean_experiment_directory()
        base = organizeData.BaseMeasurer(data_dir, work_dir, human_dir)
        print("Created checker and base for "+str(data_dir)+"; enumerating subdirectories")
        for animal_subdir in base.worm_dirs:
            if animal_subdir.split(os.path.sep)[-1] in checker.good_worms and animal_subdir.split(os.path.sep)[-1] in checker.dead_worms:
                animals_toprocess.append(animal_subdir)
        
    # Save the json with the worms of interest to run
    with (pathlib.Path(file_dir) / 'enumerated_animals.json').open('w') as subdir_file:
        json.dump({'animals_toprocess':animals_toprocess}, subdir_file)
    
    # Write the job file corresponding to this particular job
    jobscript_template = string.Template(
        '''
        #PBS -N $job_name
        #PBS -t 1-$num_toprocess
        #PBS -l nodes=1:ppn=1
        #PBS -l walltime=23:00:00
        #PBS -M plutzer@wustl.edu
        #PBS -d /scratch/plutzer/
        #PBS -e /home/plutzer/job_log_files/
        #PBS -o /home/plutzer/job_log_files/        
      	export ANACONDA_PATH="/act/Anaconda3-2.3.0/bin"
        export PATH=$$ANACONDA_PATH:$$PATH 
        source activate zplab_cluster
        cd $$PBS_O_WORKDIR
        python /home/plutzer/code/cluster_tools.py process $${PBS_ARRAYID}
        ''')
    code = jobscript_template.substitute(job_name='process_animals',
        num_toprocess=repr(len(animals_toprocess)))
    with (pathlib.Path(file_dir) / 'job_script.sh').open('w') as job_script:
        job_script.write(code)
    print("Wrote job script to "+ str(pathlib.Path(file_dir)/'job_script.sh'))    

    # Submit this job through subprocess and qsub
    #subprocess.Popen('qsub ' + str(pathlib.Path(file_dir) / 'job_script.sh'))

def process_animal(animal_id):
    with (pathlib.Path(file_dir) / 'enumerated_animals.json').open() as subdir_file:
        processing_data = json.load(subdir_file)
        
    worm_measurer = hm.HollyMeasurer(processing_data['animals_toprocess'][animal_id-1], work_dir, human_dir)
    worm_measurer.measure_a_worm()

if __name__ == "__main__":
    if sys.argv[1] == 'configure':
        configure_pipeline(*sys.argv[2:])
    elif sys.argv[1] == 'process':
        process_animal(int(sys.argv[2]))
