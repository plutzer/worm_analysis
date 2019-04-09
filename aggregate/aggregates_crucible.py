import os
import aggregates_algorithms
from skimage.io import imread
import csv


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
			gfp_filename = os.fsdecode(image)
			worm_num_ext = gfp_filename[10:]
			mask_filename = 'worm_mask_' + worm_num_ext
			gfp_image = imread(gfp_directory + '/' + gfp_filename,as_grey = True)
			mask = imread(mask_directory + '/' + mask_filename,as_grey = True)
			results = []
			results.append([worm_num_ext])
			print('Starting measurements')
			for threshold in adap_thresholds_mask:
				threshold_regions,threshold_intensity,percent_threshold = aggregates_algorithms.adap_threshold(gfp_image,threshold,mask = mask)
				results.append((threshold_regions,threshold_intensity,percent_threshold))
			print('Adap thresholds mask finsihed')
			for threshold in adap_thresholds_non:
				threshold_regions,threshold_intensity,nan = aggregates_algorithms.adap_threshold(gfp_image,threshold)
				results.append((threshold_regions,threshold_intensity))
			print('Adap thresholds non finished')
			for low in low_set_watershed_non:
				for high in high_set_watershed_non:
					results.append(aggregates_algorithms.watershed_aggs(gfp_image,low,high))
			print('Watershed non finished')
			for low in low_set_watershed_mask:
				for high in high_set_watershed_mask:
					results.append(aggregates_algorithms.watershed_aggs(gfp_image,low,high,mask))
			print('Watershed mask finished')
			for threshold in thresholds_log_dog:
				for min_sig in min_sigma:
					for max_sig in max_sigma:
						results.append(aggregates_algorithms.log(gfp_image,threshold,min_sig,max_sig))
						results.append(aggregates_algorithms.doh(gfp_image,threshold/10,min_sig,max_sig))
			print('Log, Doh finished')
			for min_sig in min_sigma:
				for max_sig in max_sigma:
					for ratio in sigma_ratio:
						for threshold in thresholds_log_dog:
							results.append(aggregates_algorithms.dog(gfp_image,min_sig,max_sig,ratio,threshold))
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

directory = '/scratch/plutzer/Aggregates_training_data/Worm1/' ## Change this as needed
run_dir(directory)

