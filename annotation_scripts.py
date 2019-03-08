from ris_widget import ris_widget; 
import elegant.gui
from elegant.gui import points_annotation
from elegant.gui import experiment_annotator
from elegant import load_data
import os
from skimage.io import imread
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