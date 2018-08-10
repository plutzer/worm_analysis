import os
import shutil


def organize_directory(worm_directory,worm_work_directory,training_directory):
	for image in os.listdir(worm_directory):
		img_filename = os.fsdecode(image)
		#print(img_filename[16:len(img_filename)-4])
		if img_filename[16:len(img_filename)-4] == 'gfp':
			timestamp = img_filename[0:15]
			shutil.copy(worm_directory + '/' + img_filename,training_directory + '/GFP/cyan_worm_' + timestamp + '.png')
			mask = worm_work_directory + '/' + timestamp + ' mask.png'
			shutil.copy(mask,training_directory + '/masks/worm_mask_' + timestamp + '.png')

