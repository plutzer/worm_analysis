import csv
from skimage.io import imread
from matplotlib import pyplot as plt
import numpy as np
from skimage.filters import sobel
from skimage.morphology import watershed
from skimage.measure import label
from scipy import stats
from skimage.feature import blob_log, blob_dog, blob_doh




def compare_annotations(annotations1,annotations2):
	annotation1_dict = {}
	annotation1 = []
	annotation2 = []
	with open(annotations1) as file:
		reader = csv.reader(file,delimiter = '\t')
		for line in reader:
			annotation1_dict[line[0]] = int(line[1])
	with open(annotations2) as file:
		reader = csv.reader(file,delimiter = '\t')
		for line in reader:
			try:
				annotation1.append(annotation1_dict[line[0]])
				annotation2.append(int(line[1]))
			except:
				do = 'nothing'
	plt.scatter(annotation1,annotation2,color='black')
	s,i,r,p,e = stats.linregress(annotation1,annotation2)
	line = []
	[line.append(a*s + i) for a in annotation1]
	plt.plot(annotation1,line,color='blue')
	top1 = np.percentile(annotation1,90)
	bottom1 = np.percentile(annotation1,10)
	top2 = np.percentile(annotation2,90)
	bottom2 = np.percentile(annotation2,10)
	plt.axvline(x=top1)
	plt.axvline(x=bottom1)
	plt.axhline(y=top2)
	plt.axhline(y=bottom2)
	[print(a) for a in [s,i,r,p,e]]
	plt.show()
	return

class tester:
	def __init__(self,annotations,directory):

		agg_function = self.blob_detection

		self.real_aggs = []
		self.measured_aggs = []
		with open(annotations) as file:
			reader = csv.reader(file,delimiter = '\t')
			for line in reader:
				self.real_aggs.append(line[1])
				image_name = line[0]
				image = directory + '/' + image_name
				im = imread(image,as_grey = True)
				self.measured_aggs.append(agg_function(im))
				print('Measured Worm:  ' + str(image_name) + '  Measured Aggregates: ' + str(self.measured_aggs[len(self.measured_aggs)-1]) + '  Annotated Aggregates: ' + str(self.real_aggs[len(self.real_aggs)-1]))
			print('--------------DONE----------------')
		self.real_aggs = np.asarray(self.real_aggs)
		self.real_aggs = self.real_aggs.astype(int)
		self.measured_aggs = np.asarray(self.measured_aggs)

	def filler(image):
		return 4

	def get_results(self):
		#ignore zeros
		self.fixed_measured_aggs = []
		self.fixed_real_aggs = []
		for item in range(0,len(self.real_aggs)):
			if self.real_aggs[item] != 0:
				self.fixed_measured_aggs.append(self.measured_aggs[item])
				self.fixed_real_aggs.append(self.real_aggs[item])
		plt.scatter(self.fixed_measured_aggs,self.fixed_real_aggs,color='black')
		slope, intercept, r_value, p_value, std_err = stats.linregress(self.fixed_measured_aggs,self.fixed_real_aggs)
		self.line = []
		[self.line.append(a*slope + intercept) for a in self.fixed_measured_aggs]
		self.line = np.asarray(self.line)
		plt.plot(self.fixed_measured_aggs,self.line,color='blue')
		plt.show()
		print([item for item in [slope,intercept,r_value, p_value, std_err]])
		self.measured_top = np.percentile(self.fixed_measured_aggs,90)
		self.measured_bottom = np.percentile(self.fixed_measured_aggs,10)
		self.real_top = np.percentile(fixed_real_aggs,90)
		self.real_bottom = np.percentile(fixed_real_aggs,10)
		counter = 0
		total = 0
		for item in range(0,len(self.fixed_measured_aggs)):
			if self.fixed_measured_aggs[item] 
				if self.fixed_measured_aggs[item] > self.measured_top and self.fixed_real_aggs[item] > self.real_top:
					counter = counter + 1
				elif self.self.fixed_measured_aggs < self.measured_bottom and self.fixed_real_aggs < self.real_bottom:
					counter = counter + 1
				total = total + 1
		same 
	def watershed_aggs(self,image):
		#Image is a GFP image as an ndarray
		high_thresh = np.percentile(image,99.99) #Alternatively, can use the 97th percentile mask pixel if the mask is available
		low_thresh = np.percentile(image,99)
		markers = np.zeros_like(image)
		markers[image > high_thresh] = 2
		markers[image < low_thresh] = 1
		watershed_im = watershed(sobel(image),markers)
		aggregates = np.max(label(watershed_im))
		return aggregates

	def integrated_percent(self,image,threshold = 99.9):
		#im = imread(image,as_grey=True)
		cutoff = np.percentile(image,threshold)
		cutoff_pix = image[image>cutoff]
		highsum = np.sum(cutoff_pix)
		totalsum = np.sum(image)
		#plt.imshow(im>cutoff)
		#plt.show()
		#print(highsum)
		#print(totalsum)
		return 100*highsum/totalsum
	
	def integrated_watershed_percent(self,image):
		high_thresh = np.percentile(image,99.99)
		low_thresh = np.percentile(image,99)
		markers = np.zeros_like(image)
		markers[image > high_thresh] = 2
		markers[image < low_thresh] = 1
		watershed_im = watershed(sobel(image),markers)
		highsum = np.sum(image[watershed_im > 1.5])
		totalsum = np.sum(image)
		return 100*highsum/totalsum

	def blob_detection(self,image):
		blobs = blob_doh(image,min_sigma = 1, max_sigma = 10, num_sigma = 20, threshold = .002, overlap = .8)
		return len(blobs)





