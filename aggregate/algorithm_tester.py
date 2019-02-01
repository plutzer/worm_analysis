import os
from skimage.io import imread
from zplab.worm_analysis.aggregate import aggregates_algorithms
import numpy as np
from matplotlib import pyplot as plt

# Get small set of test images and load them in.
def run():
	prototype_dir = os.getcwd() + '/zplab/worm_analysis' + '/testdata/prototype_data'
	for img in os.listdir(prototype_dir):
		image = imread(prototype_dir + '/' + img,as_grey = True)
		out1, out2, out3 = aggregates_algorithms.adap_threshold(img,0.15) ## Change this call as needed.
		print(args)
	return