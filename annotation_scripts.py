from ris_widget import ris_widget; 
import elegant.gui
from elegant.gui import points_annotation
from elegant.gui import experiment_annotator
from elegant.gui import pose_annotation
from elegant import load_data
from elegant import worm_spline
import os
from skimage.io import imread,imsave
import csv

## If experiment is false, returns a list of filenames and the annotations list.
## If experiment is true, returns just the annotations.
def points_annotator(directory, points_name, experiment = False): 
	rw = ris_widget.RisWidget()
	pa = points_annotation.PointsAnnotation(rw,name = points_name)
	if experiment:
		positions = load_data.scan_experiment_dir(directory)
		ea = experiment_annotator.ExperimentAnnotator(rw,directory,positions,[pa])
		return(ea)
	else:
		filenames = []
		rw.add_annotator([pa])
		for image_file in os.listdir(directory):
			image_filename = os.fsdecode(image_file)
			try:
				image = imread(directory + "/" + image_filename,as_grey = True)
				rw.flipbook_pages.append(image)
				filenames.append(image_filename)
			except:
				a = 1
		return(filenames, rw.annotator.all_annotations)

def save_annotations(annotations,filenames = [], directory = " ", experiment = False):
	if experiment:
		annotations.save_annotations()
		return
	else:
		if directory == ' ':
			print('ERROR: Need to input a directory if the annotations are not from an experiment.')
			print('Format:     save_annotations(annotations, filenames, directory, experiment)')
			return
		if len(filenames) != len(annotations):
			print('ERROR: filenames and annotations are of different lengths.')
			print('You can ignore the filenames by not passing a list of filenames to the function')
			return
		if filenames == []:
			print('Using blank filenames.')
			filenames.append(" " for item in annotations)
		with open(directory + "/" + 'annotations.txt','w') as csvfile:
			writer = csv.writer(csvfile, delimiter = '\t')
			for item in range(len(annotations)):
				writer.writerow([filenames[item],annotations[item]])
		return

def manual_straightening(directory): #Base directory: should include 'bf' and 'gfp' folders
	rw = ris_widget.RisWidget()
	pa = pose_annotation.PoseAnnotation(rw)
	rw.add_annotator([pa])
	filenames = []
	for item in os.listdir(directory + '/bf'):
		try:
			filename = os.fsdecode(item)
			fname_base = filename[0:len(filename)-6]
			#print(fname_base)
			image = imread(directory + '/bf/' + filename,as_grey = True)
			rw.flipbook_pages.append(image)
			filenames.append(filename)
		except:
			print('doing nothing')
	return(filenames,rw.annotator.all_annotations)

def create_straightened_images(directory,filenames,annotations): #Base directory: should include 'bf' and 'gfp' folders as well as 'straightened_bf' and 'straightened_gfp' folders
	for num in range(len(filenames)):
		base = filenames[num][0:len(filenames[num]) - 6]
		gfp_im = imread(directory + '/gfp/' + base + 'gfp.png',as_grey = True)
		bf_im = imread(directory + '/bf/' + base + 'bf.png',as_grey = True)
		str_bf = worm_spline.to_worm_frame(bf_im,annotations[num]['pose'][0],sample_distance = 80)
		str_gfp = worm_spline.to_worm_frame(gfp_im,annotations[num]['pose'][0],sample_distance = 80)
		imsave(directory + '/straightened_gfp/' + base + 'straightened_gfp.png',str_gfp)
		imsave(directory + '/straightened_bf/' + base + 'straightened_bf.png',str_bf)


