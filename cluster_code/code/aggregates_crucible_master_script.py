import os
from skimage.io import imread
import csv
import numpy as np
from skimage.filters import sobel
from skimage.morphology import watershed
from skimage.measure import label
from skimage.feature import blob_log, blob_dog, blob_doh

# The purpose of this script is to exhaustively test all possible aggregate counting measurements.
# Image directories with images and their associated masks will be used to make various measurements
# All of the measurements for each image will be output to a single TSV file for later analysis

'''
	Directoires must be organized in the following way:
		/directory
			/GFP
				cyan_worm_#.png
			/masks
				worm_mask_#.png
'''

def run_dir(directory):
	gfp_directory = directory + '/gfp'
	mask_directory = directory + '/masks'
	adap_thresholds_mask = [90,95,96,97,98,99,99.9,99.95,99.99]
	adap_thresholds_non = [95,97,99,99.5,99.9,99.95,99.99,100]
	low_set_watershed_non = [97,99,99.9]
	low_set_watershed_mask = [90,95,97,99]
	high_set_watershed_non = [99.9,99.95,99.99]
	high_set_watershed_mask = [97,99,99.5,99.9,99.95,99.99]
	min_sigma = [1]
	max_sigma = [4,6]
	thresholds_log_dog = [.05,.1,.125,.15,.175,.2,.25,.3,.4,.5,.6]
	thresholds_doh = [item/10 for item in thresholds_log_dog]
	sigma_ratio = [1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,5,10]

	with open(directory + '/aggregate_measurements_save.tsv','w') as csvfile:
		measurement_writer = csv.writer(csvfile,delimiter = '\t')
		headers = []
		headers.append('worm_num')
		for threshold in adap_thresholds_mask:
			for measurement in ['t_reg_mask_','t_int_mask_','t_perc_mask_']:
				headers.append(measurement + str(threshold))
		for threshold in adap_thresholds_non:
			for measurement in ['t_reg_non_','t_int_non_']:
				headers.append(measurement + str(threshold))
		for low in low_set_watershed_non:
			for high in high_set_watershed_non:
				for measurement in ['wat_reg_non_','wat_per_non_']:
					headers.append(measurement + str(low) + '_' + str(high))
		for low in low_set_watershed_mask:
			for high in high_set_watershed_mask:
				for measurement in ['wat_reg_mask_','wat_per_mask_']:
					headers.append(measurement + str(low) + '_' + str(high))
		for threshold in thresholds_log_dog:
			for max_sig in max_sigma:
				headers.append('log_' + str(threshold) + '_' + str(max_sig))
				headers.append('doh_' + str(threshold/10) + '_' + str(max_sig))
		for max_sig in max_sigma:
			for ratio in sigma_ratio:
				for threshold in thresholds_log_dog:
					headers.append('dog_' + str(max_sig) + '_' + str(ratio) + '_' + str(threshold))
		measurement_writer.writerow([item for item in headers]) #Change this	
		for image in os.listdir(gfp_directory):
			try:
				gfp_filename = os.fsdecode(image)
				worm_num_ext = gfp_filename[10:]
				mask_filename = 'worm_mask_' + worm_num_ext
				gfp_image = imread(gfp_directory + '/' + gfp_filename,as_grey = True)
				mask = imread(mask_directory + '/' + mask_filename,as_grey = True)
				results = []
				results.append([worm_num_ext])
				print('Starting measurements')
				for threshold in adap_thresholds_mask:
					threshold_regions,threshold_intensity,percent_threshold = adap_threshold(gfp_image,threshold,mask = mask)
					results.append((threshold_regions,threshold_intensity,percent_threshold))
				print('Adap thresholds mask finsihed')
				for threshold in adap_thresholds_non:
					threshold_regions,threshold_intensity,nan = adap_threshold(gfp_image,threshold)
					results.append((threshold_regions,threshold_intensity))
				print('Adap thresholds non finished')
				for low in low_set_watershed_non:
					for high in high_set_watershed_non:
						results.append(watershed_aggs(gfp_image,low,high))
				print('Watershed non finished')
				for low in low_set_watershed_mask:
					for high in high_set_watershed_mask:
						results.append(watershed_aggs(gfp_image,low,high,mask))
				print('Watershed mask finished')
				for threshold in thresholds_log_dog:
					for min_sig in min_sigma:
						for max_sig in max_sigma:
							results.append(log(gfp_image,threshold,min_sig,max_sig))
							results.append(doh(gfp_image,threshold/10,min_sig,max_sig))
				print('Log, Doh finished')
				for min_sig in min_sigma:
					for max_sig in max_sigma:
						for ratio in sigma_ratio:
							for threshold in thresholds_log_dog:
								results.append(dog(gfp_image,min_sig,max_sig,ratio,threshold))
				items_list = []
				print('Dog finished')
				for measurement in results:
					try:
						for item in measurement:
							items_list.append(item)
					except:
						items_list.append(measurement)
				measurement_writer.writerow([item for item in items_list])
				print('Wrote results for image: ' + str(worm_num_ext))
			except:
				print('File analysis failed:   ' + str(gfp_filename))

def adap_threshold(image,thresh,mask = None):
	if mask is None:
		threshold = np.percentile(image,thresh)
	else:
		inMask = image[mask > .5]
		threshold = np.percentile(inMask,thresh)
	threshIm = image>threshold
	labeled_threshIm = label(threshIm)
	threshold_regions = np.max(labeled_threshIm)
	threshold_intensity = threshold
	if mask is not None:
		percent_threshold = (np.sum(threshIm))/(np.sum(mask > .5))
		return(threshold_regions,threshold_intensity,percent_threshold)
	else:
		return(threshold_regions,threshold_intensity,'NaN')

def watershed_aggs(image,low,high,mask = None):
	if mask is None:
		high_thresh = np.percentile(image,high)
		low_thresh = np.percentile(image,low)
	else:
		high_thresh = np.percentile(image[mask>.5],high)
		low_thresh = np.percentile(image[mask>.5],low)
	markers = np.zeros_like(image)
	markers[image>high_thresh] = 2
	markers[image<low_thresh] = 1
	watershed_im = watershed(sobel(image),markers)
	watershed_regions = np.max(label(watershed_im))-1
	highsum = np.sum(image[watershed_im > 1.5])
	totalsum = np.sum(image)
	integrated_watershed_percent = 100*highsum/totalsum
	return(watershed_regions,integrated_watershed_percent)


def log(image,thresh,min_sigma,max_sigma,num_sigma = 8,mask = None):
	blobs = blob_log(image,min_sigma = min_sigma,max_sigma = max_sigma,num_sigma = num_sigma,threshold = thresh,overlap = .8)
	total_blobs = len(blobs)
	return(total_blobs)

def doh(image,thresh,min_sigma,max_sigma,num_sigma = 8,mask = None):
	blobs = blob_doh(image,min_sigma = min_sigma,max_sigma = max_sigma,num_sigma = num_sigma,threshold = thresh,overlap = .8)
	total_blobs = len(blobs)
	return(total_blobs)

def dog(image,min_sigma,max_sigma,sigma_ratio,thresh,mask = None):
	blobs = blob_dog(image,min_sigma = min_sigma,max_sigma = max_sigma,sigma_ratio = sigma_ratio,threshold = thresh,overlap = .8)
	total_blobs = len(blobs)
	return(total_blobs)


directories = ['/scratch/plutzer/Aggregates_training_data/sorter_images','/scratch/plutzer/Aggregates_training_data/Worm1','/scratch/plutzer/Aggregates_training_data/Worm2','/scratch/plutzer/Aggregates_training_data/Worm3','/scratch/plutzer/Aggregates_training_data/Worm4']
for directory in directories:
        print('Running directory:   ' + str(directory))
        run_dir(directory)


