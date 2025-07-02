import matplotlib.pyplot as plt  # standard Python plotting library
import numpy as np  # fundamental package for scientific computing, handles arrays and math
import tdt # import the tdt library
import os
import pandas as pd

file_paths = [
    #(r"C:\Users\stuberadmin\Downloads\red code\0007-0002-250401-113301", 624),
    (r"C:\Users\stuberadmin\Downloads\red code\0007-0001-250401-103502")

]

channels_to_include = ['redB', 'isoA']

### don't edit below

df_column_names = ['', 'blockname', 'start_date', 'channel', 'raw_au']
block_path = file_paths[0]

# Read all stream data into a Python structure called 'dataft'
dataft = tdt.read_block(block_path, evtype=['streams'])
block_name = dataft['info']['blockname']

save_path = os.path.join(block_path, f'{block_name}_streams_data')

# Filter out channels that do not start with 'Fi'
channels = [k for k, v in dataft.streams.items() if not k.startswith('Fi')]
if channels_to_include:
    channels = list(set(channels_to_include) & set(channels))

df_data = pd.DataFrame(columns = list(range(len(df_column_names))))
df_data = df_data.reset_index(drop=True)
for channel in channels:
    
    channel_data = dataft['streams'][channel]['data']
    
    df_tmp = pd.DataFrame({
            0: list(range(len(channel_data))),
            1: [block_name] * len(channel_data),
            2: [dataft['info']['start_date'] ] * len(channel_data),
            3: [channel] * len(channel_data),
            4: channel_data
        })
    
    df_data = pd.concat([df_data, df_tmp], axis=0)

#df_data = df_data.reset_index(drop=True)
df_data = df_data.astype({0: 'int32', 1: 'string', 2: 'string', 3: 'string', 4: 'float32'})      
 

#set column names equal to values in row index position 0
df_data.columns = df_column_names

df_data.to_feather(save_path + '.feather')
df_data.to_csv(save_path + '.csv', header=False, index=False)

# Process each file and store the results
# for file_path in file_paths:
#     processed_channels = process_file(file_path)  # Change to 'linear' if needed