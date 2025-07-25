from re import search
from tdt import read_block
import pandas as pd
import os
import shutil
import tempfile

# info -----------------------------------------------------------------------------------------------------------------
def tidy_tdt_info(data_tdt):
    # extracts info from tdt structure

    tidy_info = pd.DataFrame(
        {
            'blockname': [data_tdt.info.blockname],
            'tankpath': [data_tdt.info.tankpath],
            'start_date': [data_tdt.info.start_date],
            'utc_start_time': [data_tdt.info.utc_start_time],
            'stop_date': [data_tdt.info.stop_date],
            'utc_stop_time': [data_tdt.info.utc_stop_time],
            'duration': [data_tdt.info.duration],
            'stream_channel': [data_tdt.info.stream_channel],
            'snip_channel': [data_tdt.info.snip_channel]

        })

    return tidy_info

# streams --------------------------------------------------------------------------------------------------------------
def extract_stream_info(data_tdt_streams):
    # return info for an individual tdt stream

    tidy_streams_info = pd.DataFrame(
        {
            'name': data_tdt_streams.name,
            'code': data_tdt_streams.code,
            'size': [data_tdt_streams.size],
            'type': data_tdt_streams.type,
            'type_str': data_tdt_streams.type_str,
            'ucf': data_tdt_streams.ucf,
            'fs': data_tdt_streams.fs,
            'dform': data_tdt_streams.dform,
            'start_time': data_tdt_streams.start_time,
            'channel': data_tdt_streams.channel
        })

    return (tidy_streams_info)


def extract_stream_data(data_tdt_streams):
    # return data for an individual tdt stream

    tidy_streams_data = pd.DataFrame({
        'channel': data_tdt_streams.name,
        'raw_au': data_tdt_streams.data
    })

    return (tidy_streams_data)


def tidy_tdt_streams(data, channel_names):
    # returns data and info from all streams within tdt structure

    first_concat = 1

    for stream in data.streams.keys():
        if (search('_\d\d\d[A-Za-z]', stream) or stream in channel_names):  # if stream contains _###C (id for fp channels; filters tdt preprocessed streams)
            if (first_concat):
                streams_info = extract_stream_info(data.streams[stream])
                streams_data = extract_stream_data(data.streams[stream])

                first_concat = 0;
            else:
                streams_info = pd.concat([streams_info, extract_stream_info(data.streams[stream])])
                streams_data = pd.concat([streams_data, extract_stream_data(data.streams[stream])])

    streams_info.insert(0, 'start_date', data.info.start_date)
    streams_info.insert(0, 'blockname', data.info.blockname)

    streams_data.insert(0, 'start_date', data.info.start_date)
    streams_data.insert(0, 'blockname', data.info.blockname)

    return (streams_info, streams_data)

# epochs ---------------------------------------------------------------------------------------------------------------
def extract_epoch_info(data_tdt_streams):
    # return info for an individual tdt epoc

    tidy_epoch_info = pd.DataFrame(
        {
            'name': data_tdt_streams.name,
            'type': data_tdt_streams.type,
            'type_str': data_tdt_streams.type_str,
            'dform': data_tdt_streams.dform,
            'size': [data_tdt_streams.size]
        })

    return (tidy_epoch_info)


def extract_epoch_data(data_tdt_streams):
    # return data for an individual tdt epoc

    tidy_epoch_data = pd.DataFrame(
        {'name': data_tdt_streams.name,
         'onset': data_tdt_streams.onset,
         'offset': data_tdt_streams.offset
         })

    return (tidy_epoch_data)


def tidy_tdt_epocs(data):
    flag_epoch = 0

    for epoch in data.epocs.keys():
        flag_epoch = 1

        if (epoch == list(data.epocs.keys())[0]):
            epocs_info = extract_epoch_info(data.epocs[epoch])
            epocs_data = extract_epoch_data(data.epocs[epoch])

        else:
            epocs_info = pd.concat([epocs_info, extract_epoch_info(data.epocs[epoch])])
            epocs_data = pd.concat([epocs_data, extract_epoch_data(data.epocs[epoch])])

    if flag_epoch:
        epocs_info.insert(0, 'start_date', data.info.start_date)
        epocs_info.insert(0, 'blockname', data.info.blockname)

        epocs_data.insert(0, 'start_date', data.info.start_date)
        epocs_data.insert(0, 'blockname', data.info.blockname)

        return (epocs_info, epocs_data, flag_epoch)
    else:
        return(0,0, flag_epoch)

# extract and tidy all files in raw directory that are not located in extracted directory ------------------------
def tidy_tdt_extract_and_tidy(dir_raw, dir_extracted, channel_names):
    # return lists of files in raw and processed directories___
    raw_block_paths = os.listdir(dir_raw)
    processed_files = os.listdir(dir_extracted)

    # filter out hidden files from raw_block_paths
    raw_block_paths = list(filter(lambda raw_block_paths: raw_block_paths.find('.') == -1, raw_block_paths))

    # trim file names to blocknames
    processed_files = {x.replace('_streams_info.feather','') for x in processed_files}
    processed_files = {x.replace('_streams_data.feather','') for x in processed_files}
    processed_files = {x.replace('_epocs_info.feather','')   for x in processed_files}
    processed_files = {x.replace('_epocs_data.feather','')   for x in processed_files}
    processed_files = {x.replace('_info.feather','')         for x in processed_files}

    processed_files = list(processed_files)

    # remove files that have already been processed____________
    process_block_paths = raw_block_paths

    for processed_file in processed_files:
        try:
            while True:
                process_block_paths.remove(processed_file)
        except ValueError:
            pass


    if len(process_block_paths) >= 1:
        print('')
        print('extracting data from dir: ' + dir_raw + ' to dir: ' + dir_extracted)

        for process_block_path in process_block_paths:
            block_path = os.path.join(dir_raw, process_block_path)

            print("extracting blockpath: " + block_path)

            data = read_block(block_path, evtype = ['all'])

            data_info = tidy_tdt_info(data)
            streams_info, streams_data = tidy_tdt_streams(data, channel_names)
            epocs_info, epocs_data, flag_epoch = tidy_tdt_epocs(data)

            session_id = data_info['blockname'][0]

            data_info.to_feather(os.path.join(dir_extracted, session_id +'_info.feather'))

            streams_info.reset_index(drop = True).to_feather(os.path.join(dir_extracted, session_id + '_streams_info.feather'))
            streams_data.reset_index(drop = True).to_feather(os.path.join(dir_extracted, session_id + '_streams_data.feather'))

            if(flag_epoch):
                epocs_info.reset_index(drop = True).to_feather(os.path.join(dir_extracted, session_id + '_epocs_info.feather'))
                epocs_data.reset_index(drop = True).to_feather(os.path.join(dir_extracted, session_id + '_epocs_data.feather'))
                
            print("Done with blockpath: " + block_path)
    else:
        print('no files to extract... all fp in dir :'+ dir_raw + ' has already been extracted to dir: ' + dir_extracted)

#doric ---------------------------------------------------------------------------------------------------------
def doric_convert_epoch_id(name):
    if name == 'DIO01':
        name = 'PtC0'

    if name == 'DIO02':
        name = 'PtC1'

    if name == 'DIO03':
        name = 'PtC2'

    if name == 'DIO04':
        name = 'PtC3'

    return (name)

def doric_extract_epoch_data(data):
    import numpy as np

    for stream in data:
        name = stream["Name"]

        # filter data to DigitalIO
        if 'DigitalIO' in stream["Name"]:
            stream_digital_io = stream

            # save time
            for io in stream_digital_io["Data"]:
                if 'Time' in io["Name"]:
                    io_time = io["Data"]

            first_concat = 1

            # save onset times and save as a data frame
            for io in stream_digital_io["Data"]:
                if 'Time' not in io["Name"]:
                    name = io["Name"]
                    name = doric_convert_epoch_id(name)


                    io_state = io["Data"]

                    idx_onset = np.where(np.diff(io_state) == 1)[0] + 1
                    onset = io_time[idx_onset]

                    idx_offset = np.where(np.diff(io_state) == -1)[0]
                    offset = io_time[idx_offset]

                    io_df = pd.DataFrame(
                        {
                            'name': name,
                            'onset': onset,
                            'offset': offset
                        })

                    if first_concat:
                        io_df_combined = io_df
                        first_concat = 0

                    else:
                        io_df_combined = pd.concat([io_df_combined, io_df])

    return (io_df_combined)

def doric_convert_channel_id(channel_id):
    if channel_id == 'AIN01xAOUT01':
        channel_id = '_560B'

    if channel_id == 'AIN02xAOUT01':
        channel_id = '_470A'

    if channel_id == 'AIN02xAOUT02':
        channel_id = '_405A'

    if channel_id == 'AIN03xAOUT01':
        channel_id = '_560D'

    if channel_id == 'AIN04xAOUT01':
        channel_id = '_470C'

    if channel_id == 'AIN04xAOUT02':
        channel_id = '_405C'

    return (channel_id)


def doric_extract_stream_info(data):
    # return info for an individual doric stream
    import numpy as np
    
    streams_info_list = []
    
    for stream in data:
        name = stream["Name"]
        
        # filter data to LockIn streams
        if 'LockIn' in stream["Name"]:
            channel_id = name.split('_')[-1]
            channel_id = channel_id.replace('-LockIn', '')
            channel_id = doric_convert_channel_id(channel_id)
            
            # Extract time data to calculate sampling rate
            time_data = None
            for io in stream["Data"]:
                if 'Time' in io["Name"]:
                    time_data = io["Data"]
                    break
            
            # Calculate sampling rate if time data is available
            fs = None
            if time_data is not None and len(time_data) > 1:
                # Calculate average sampling rate from time differences
                time_diffs = np.diff(time_data)
                fs = 1.0 / np.mean(time_diffs) if np.mean(time_diffs) > 0 else None
            
            # Get data size
            data_size = len(time_data) if time_data is not None else 0
            
            # Create stream info DataFrame
            stream_info = pd.DataFrame({
                'name': [channel_id],
                'channel': 1,
                'fs': [fs],
                'size': [data_size],
                'type': ['LockIn'],
                'type_str': ['LockIn Stream'],
                'start_time': [time_data[0] if time_data is not None and len(time_data) > 0 else None],
                'end_time': [time_data[-1] if time_data is not None and len(time_data) > 0 else None]
            })
            
            streams_info_list.append(stream_info)
    
    if streams_info_list:
        return pd.concat(streams_info_list, ignore_index=True)
    else:
        return pd.DataFrame()


def doric_extract_stream_data(data):
    first_concat = 1

    for stream in data:
        name = stream["Name"]

        # filter data to LockIn
        if 'LockIn' in stream["Name"]:

            channel_id = name.split('_')[-1]
            channel_id = channel_id.replace('-LockIn', '')
            channel_id = doric_convert_channel_id(channel_id)

            # save time
            for io in stream["Data"]:
                if 'Time' in io["Name"]:
                    time = io["Data"]

            # save channel_id

            # save onset times and save as a data frame
            for io in stream["Data"]:
                if 'Time' not in io["Name"]:
                    raw_au = io["Data"]

                    stream_df = pd.DataFrame(
                        {
                            'time': time,
                            'channel': channel_id,
                            'raw_au': raw_au
                        })

                    if first_concat:
                        stream_df_combined = stream_df
                        first_concat = 0

                    else:
                        stream_df_combined = pd.concat([stream_df_combined, stream_df])

    return (stream_df_combined)

def tidy_doric_extract_and_tidy(dir_raw, dir_extracted):
    import doric as dr

    # return lists of files in raw and processed directories___
    raw_block_paths = os.listdir(dir_raw)
    processed_files = os.listdir(dir_extracted)

    # trim file names to blocknames
    raw_block_paths = {x.replace('.doric','') for x in raw_block_paths}

    # trim file names to blocknames
    processed_files = {x.replace('_streams_data.feather','') for x in processed_files}
    processed_files = {x.replace('_epocs_data.feather','')   for x in processed_files}

    processed_files = list(processed_files)

    # remove files that have already been processed____________
    process_block_paths = raw_block_paths


    for processed_file in processed_files:
        try:
            while True:
                process_block_paths.remove(processed_file)
        except KeyError or ValueError:
            pass


    if len(process_block_paths) >= 1:
        print('')
        print('extracting data from dir: ' + dir_raw + ' to dir: ' + dir_extracted)

        for process_block_path in process_block_paths:
            block_path = os.path.join(dir_raw, process_block_path)

            print("extracting blockpath: " + block_path)

            session_id = process_block_path.replace('.doric', '')

            data = dr.ExtractDataAcquisition(block_path + '.doric') # read in data using doric function

            streams_info = doric_extract_stream_info(data)
            streams_info["blockname"] = session_id
            
            streams_data = doric_extract_stream_data(data)
            streams_data["blockname"] = session_id

            epocs_data = doric_extract_epoch_data(data)
            epocs_data["blockname"] = session_id

            streams_info.reset_index(drop = True).to_feather(os.path.join(dir_extracted, session_id + '_streams_info.feather'))
            streams_info.reset_index(drop = True).to_csv(os.path.join(dir_extracted, session_id + '_streams_info.csv'))
            streams_data.reset_index(drop = True).to_feather(os.path.join(dir_extracted, session_id + '_streams_data.feather'))
            streams_data.reset_index(drop = True).to_csv(os.path.join(dir_extracted, session_id + '_streams_data.csv'))
            epocs_data.reset_index(drop = True).to_feather(os.path.join(dir_extracted, session_id + '_epocs_data.feather'))
    else:
        print('no files to extract... all fp in dir :'+ dir_raw + ' has already been extracted to dir: ' + dir_extracted)
        
        
if __name__ == "__main__":
    dir_raw = r'D:\photom\raw'
    dir_extracted = r'D:\photom\extracted'
    channel_names = ['405A', '465A'] # not needed for doric

    doric_files = []
    tdt_block_dirs = set()

    # Recursively walk through dir_raw to find .doric and .tev files
    for root, dirs, files in os.walk(dir_raw):
        for file in files:
            if file.lower().endswith('.doric'):
                doric_files.append(os.path.join(root, file))
            elif file.lower().endswith('.tev'):
                tdt_block_dirs.add(root)

    # Process doric files if any
    if doric_files:
        with tempfile.TemporaryDirectory() as temp_doric_dir:
            for f in doric_files:
                shutil.copy(f, os.path.join(temp_doric_dir, os.path.basename(f)))
            tidy_doric_extract_and_tidy(temp_doric_dir, dir_extracted)

    # Process TDT block folders if any
    if tdt_block_dirs:
        with tempfile.TemporaryDirectory() as temp_tdt_dir:
            for block_dir in tdt_block_dirs:
                dst = os.path.join(temp_tdt_dir, os.path.basename(block_dir))
                if os.path.isdir(block_dir):
                    shutil.copytree(block_dir, dst)
            tidy_tdt_extract_and_tidy(temp_tdt_dir, dir_extracted, channel_names)