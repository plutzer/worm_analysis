import pathlib
import shutil
import json
import numpy as np

import matplotlib.pyplot as plt

def remove_offending_tp(expt_path,timept_str,dry_run=False):
    '''
        expt_path - str/pathlib.Path pointing to expt. dir
        timept_str - str in format yyyymmdd-thhmm
        dry_run - Toggles taking action (if False, will not delete but will verbosely specify where offending files are found
    '''
    
    if type(expt_path) is str: expt_path = pathlib.Path(expt_path)
    
    for sub_dir in expt_path.iterdir():
        if sub_dir.is_dir():
            
            # Move bad files
            offending_files = [im_file for im_file in sub_dir.iterdir() if timept_str in str(im_file)]
            if len(offending_files)>0:
                print('Found offending files in: '+str(sub_dir))
                
                if not dry_run:
                    if not (sub_dir/'offending_files').exists(): (sub_dir/'offending_files').mkdir(mode=755)
                    [im_file.rename(im_file.parent/'offending_files'/im_file.parts[-1]) for im_file in offending_files]
            
            # Check if position metadata is bad and handle 
            if sub_dir.parts[-1] != 'calibrations':
                md_file = (sub_dir/'position_metadata.json')
                if not md_file.exists(): print('No metadata file found at:'+str(sub_dir))
                else:
                    with md_file.open() as md_fp:
                        pos_md = json.load(md_fp)
                    has_bad_md_tmpt = any([tmpt['timepoints'] == timept_str for tmpt in pos_md])
                    if has_bad_md_tmpt:
                        print('Found offending entry in position_metadata in: '+str(sub_dir))
                        if not dry_run:
                            if not (sub_dir/'offending_files').exists(): (sub_dir/'offending_files').mkdir(mode=755)
                            md_file.rename(md_file.parent/'offending_files'/'position_metadata_old.json')     # Backup old
                            pos_md = [tp_data for tp_data in pos_md if tp_data['timepoints'] != timept_str]
                            with (sub_dir/'position_metadata.json').open('w') as md_fp:
                                encode_legible_to_file(pos_md,md_fp)   # Write out new position_metadata
    
    md_file = (expt_path/'experiment_metadata.json')
    with md_file.open() as md_fp:
        expt_md = json.load(md_fp)
    
    try:
        tp_idx = np.where(np.array(expt_md['timepoints']) == timept_str)[0][0]
        del expt_md['timepoints'][tp_idx]
        del expt_md['timestamps'][tp_idx]
        del expt_md['durations'][tp_idx]
        expt_md['brightfield_metering']={key:val for key,val in expt_md['brightfield_metering'] if key != timept_str}
        
        md_file.rename(md_file.parent/'experiment_metadata_old.json') #Backup
        with (expt_path/'experiment_metadata.json').open('w') as md_fp:  # Write out new
            encode_legible_to_file(expt_md,md_fp)
    except:  # Offending timepoint didn't make it into metadata
        pass

#~ def reset_mds(expt_path):
    #~ if type(expt_path) is str: expt_path = pathlib.Path(expt_path)  
    #~ for sub_dir in expt_path.iterdir():
        #~ if sub_dir.is_dir() and (sub_dir/'position_metadata_old.json').exists():
            #~ (sub_dir/'position_metadata_old.json').rename(sub_dir/'position_metadata.json')
            #~ (sub_dir/'position_metadata_old.json').unlink()

def check_expt_movement_acquisitions(expt_dir, mode='absolute'):
    if type(expt_dir) is str: expt_dir = pathlib.Path(expt_dir)
    expt_pos_data, dirs_with_mdata, dirs_miss_mdata = load_movement_metadata(expt_dir)

    expt_mdata = json.load((expt_dir/'experiment_metadata.json').open())
    if mode == 'absolute':  # Look at the absolut time at which a frame was taken (with 0 being the first bright field image)
        unpacked_times = [[[item[i] for item in expt_pos_data['movement_frame_times'][timepoint]] for timepoint in expt_mdata['timepoints']] for i in range(5)] # unpacked_times [=] [pic_number, timepoint, well_number]

        fig_h,ax_h = plt.subplots(1,5)
        [[movement_ax.scatter(idx*np.ones([len(timepoint_times),1]),np.array(timepoint_times)) for idx,timepoint_times in enumerate(frame_times)] for movement_ax,frame_times in zip(ax_h,unpacked_times)] 
        [movement_ax.set_title('Movement Frame {}'.format(n+1)) for n,movement_ax in enumerate(ax_h)]
        [movement_ax.set_xlabel('Timepoint Index') for movement_ax in ax_h]
        [movement_ax.set_ylabel('Time Taken (s)') for movement_ax in ax_h]
    elif mode == 'deltas': # Look at time taken between each set of consecutive frames in an acquisition
        expt_pos_data_diff = {timepoint:[np.diff(item) for item in expt_pos_data['movement_frame_times'][timepoint]] for timepoint in expt_mdata['timepoints']} #{timepoint: list of lists}
        unpacked_diffs = [[[item[i] for item in expt_pos_data_diff[timepoint]] for timepoint in expt_mdata['timepoints']] for i in range(4)] # unpacked_times [=] [pic_number, timepoint, well_number]
        
        fig_h,ax_h = plt.subplots(1,4)
        [[movement_ax.scatter(idx*np.ones([len(timepoint_times),1]),np.array(timepoint_times)) for idx,timepoint_times in enumerate(frame_diffs)] for movement_ax,frame_diffs in zip(ax_h,unpacked_diffs)] 
        [movement_ax.set_title('Interval b/t Frames {} & {}'.format(n+1,n+2)) for n,movement_ax in enumerate(ax_h)]
        [movement_ax.set_xlabel('Timepoint Index') for movement_ax in ax_h]
        [movement_ax.set_ylabel('Time Taken (s)') for movement_ax in ax_h]

def copy_position_md(expt_dir, dest_dir):
    '''
        Copies position metadatas all to one destination folder ('dest_dir') for backing up
    '''
    
    if type(expt_dir) is str: expt_dir = pathlib.Path(expt_dir)
    if type(dest_dir) is str: dest_dir = pathlib.Path(dest_dir)
    for sub_dir in expt_dir.iterdir():
        if sub_dir.is_dir() and str(sub_dir.parts[-1]) != 'calibrations':
            if not (dest_dir/sub_dir.parts[-1]).exists(): (dest_dir/sub_dir.parts[-1]).mkdir(mode=744)
            shutil.copyfile(str(sub_dir/'position_metadata.json'),str((dest_dir/sub_dir.parts[-1])/'position_metadata.json'))


def impute_missing_metadata_WZacquisition(expt_dir, dry_run = False):
    if type(expt_dir) is str: expt_dir = pathlib.Path(expt_dir)
    expt_pos_data, dirs_with_mdata, dirs_miss_mdata = load_movement_metadata(expt_dir)  #(expt_pos_data = {data_key: timepoint: list of list of values})
    
    expt_mdata = json.load((expt_dir/'experiment_metadata.json').open())
    unpacked_times_tmpt = [[[item[pic_num] for item in expt_pos_data['movement_frame_times'][timepoint]] for pic_num in range(5)] for timepoint in expt_mdata['timepoints']]# unpacked_times_tmpt [=] [timepoint,pic_number,well_number]
    
    entry_dict = {
        'fine_z':None,  # Not needed
        'image_timestamps': {
            'bf.png': 0.0
        },
        'movement_frame_times': None,
        'timepoint':None,
        'timepoints':None,
        'timestamp':None,
        'notes':'added by automatic imputation (DS)',
    }
    # WARNING!!!!!! TODO Make this scalable to arbitrary position metadata....
    
    new_md_data = []
    for timepoint,timestamp,timepoint_times in zip(expt_mdata['timepoints'],expt_mdata['timestamps'],unpacked_times_tmpt):
        new_entry = entry_dict.copy()
        new_entry['timepoint'] = new_entry['timepoints'] = timepoint
        new_entry['timestamp'] = timestamp
        new_entry['movement_frame_times'] =[np.median(frame_times) for frame_times in timepoint_times]
        new_md_data.append(new_entry)
    
    if dry_run:
        print('New metadata contents:')
        print(new_md_data)
        
        print('Directories missing metadata')
        print(dirs_miss_mdata)
    else:
        for pos_dir in dirs_miss_mdata:
            print('Writing new metadata to directory: '+str(pos_dir))
            with (pos_dir/'position_metadata.json').open('w') as md_fp:
                encode_legible_to_file(new_md_data,md_fp)

def load_movement_metadata(expt_dir):
    if type(expt_dir) is str: expt_dir = pathlib.Path(expt_dir)

    # Load experiment_metadata and grab timepoints to populate
    expt_mdata = json.load((expt_dir/'experiment_metadata.json').open())


    # Find dirs with position_metadata
    dirs_with_mdata = []
    dirs_miss_mdata = []
    for pos_dir in expt_dir.iterdir():
        if (pos_dir/'position_metadata.json').exists(): dirs_with_mdata.append(pos_dir)
        elif pos_dir.is_dir() and str(pos_dir.parts[-1]) != 'calibrations': dirs_miss_mdata.append(pos_dir)
    
    pos_mdata = json.load((dirs_with_mdata[0] / 'position_metadata.json').open())
    mdata_keys = [mdata_key for mdata_key in pos_mdata[0].keys()]
    expt_pos_data = {mdata_key:{timepoint:[] for timepoint in expt_mdata['timepoints']} \
        for mdata_key in mdata_keys}    #(expt_pos_data = {data_key: timepoint: list of values})

    for pos_dir in dirs_with_mdata:
        pos_mdata = json.load((pos_dir/'position_metadata.json').open())
        for item in pos_mdata:
            for mdata_key in mdata_keys:
                expt_pos_data[mdata_key][item['timepoint']].append(item[mdata_key])
    
    return (expt_pos_data, dirs_with_mdata, dirs_miss_mdata)
    
'''
age-1 run 22
bad timepoint: 20160713-1456
position_metadata_old still in 00 (still with entry for offending timepoint)
position_metadata in 64-84; offending timepoint purged

'''


# Swiped from zplab.rpc_acquisition....json_encode
class Encoder(json.JSONEncoder):
    """JSON encoder that is smart about converting iterators and numpy arrays to
    lists, and converting numpy scalars to python scalars.
    Caution: it is absurd to send large numpy arrays over the wire this way. Use
    the transfer_ism_buffer tools to send large data.
    """
    def default(self, o):
        try:
            return super().default(o)
        except TypeError as x:
            if isinstance(o, np.generic):
                item = o.item()
                if isinstance(item, npy.generic):
                    raise x
                else:
                    return item
            try:
                return list(o)
            except:
                raise x


COMPACT_ENCODER = Encoder(separators=(',', ':'))
READABLE_ENCODER = Encoder(indent=4, sort_keys=True)

def encode_compact_to_bytes(data):
    return COMPACT_ENCODER.encode(data).encode('utf8')

def encode_legible_to_file(data, f):
    for chunk in READABLE_ENCODER.iterencode(data):
        f.write(chunk)
