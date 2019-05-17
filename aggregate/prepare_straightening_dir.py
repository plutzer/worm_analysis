import os
import shutil

def organize_dir(directory,new_directory): #Directory is the experiment directory and new_directory is a folder containting 'bf' and 'gfp' folders
	bf_dir = new_directory + '/bf'
	gfp_dir = new_directory + '/gfp'
	for item in os.listdir(directory):
		item_filename = os.fsdecode(item)
		#print(item_filename)
		if len(item_filename) > 7 :
			filetag = item_filename[len(item_filename)-7:len(item_filename)-4]
		else:
			filetag = ' '
		if filetag == 'gfp':
			fname_base = item_filename[0:len(item_filename)-7]
			brightfield_fname = fname_base + 'bf.png'
			shutil.copy(directory + '/' + item_filename, gfp_dir + '/' + item_filename)
			try:
				shutil.copy(directory + '/' + brightfield_fname, bf_dir + '/' + brightfield_fname)
			except:
				print('Could not find the accompanying brightfield image for: ' + item_filename)
				os.remove(gfp_dir + '/' + item_filename)



