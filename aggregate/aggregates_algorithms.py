from skimage.io import imread
from matplotlib import pyplot as plt
import numpy as np
from skimage.filters import sobel
from skimage.morphology import watershed
from skimage.measure import label
from skimage.feature import blob_log,blob_dog,blob_doh

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

