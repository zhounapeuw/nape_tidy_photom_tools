import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

def load_data(basepath):
    for ext in ['.feather', '.csv']:  # Prefer feather if both exist
        file = basepath + ext
        if os.path.exists(file):
            if ext == '.feather':
                return pd.read_feather(file)
            else:
                return pd.read_csv(file)
    raise FileNotFoundError(f"No CSV or Feather file found for base name: {basepath}")

def fp_moving_average(streams_data, streams_info, smoothing_window) -> pd.DataFrame:
    """
    Returns moving average of raw_au as au_smooth_ma.

    Parameters:
    - streams_data: DataFrame from *_streams_data.feather exported from Python extractor
      (fp data stream variable name needs to be raw_au).
    - streams_info: DataFrame from *_streams_info.feather exported from Python extractor.
    - smoothing_window: Duration of window (seconds) for moving average to be applied.
      A smoothing window of 0 or NA will return streams_data without applying a moving average.

    Returns:
    - DataFrame with smoothed variable if smoothing_window > 0, else returns the original DataFrame.
    """
    if smoothing_window <= 0 or pd.isna(smoothing_window):
        return streams_data  # Return original DataFrame without smoothing

    # Extract channel info including sampling rate
    channel_info = streams_info[['blockname', 'name', 'fs']].drop_duplicates()

    streams_data_smoothed = pd.DataFrame()

    # Calculate moving average for each channel independently
    for _, row in channel_info.iterrows():
        blockname = row['blockname']
        channel = row['name']
        fs = pd.to_numeric(row['fs'], errors='coerce')

        if pd.isna(fs) or fs <= 0:
            continue  # Skip invalid sampling rates

        # Filter streams_data to a single channel
        streams_data_loop = streams_data[
            (streams_data['blockname'] == blockname) &
            (streams_data['channel'] == channel)
        ].copy()

        # Calculate smoothing window in samples
        smoothing_window_samples = round(fs * smoothing_window)

        # Apply centered moving average with NA padding
        streams_data_loop['au_smooth_ma'] = (
            streams_data_loop['raw_au']
            .rolling(window=smoothing_window_samples, center=True, min_periods=1)
            .mean()
        )

        # Append to result
        streams_data_smoothed = pd.concat([streams_data_smoothed, streams_data_loop], ignore_index=True)

    # Filter out rows where smoothing could not be computed (to match R behavior)
    streams_data_smoothed = streams_data_smoothed.dropna(subset=['au_smooth_ma'])

    return streams_data_smoothed

fs = 1017.25262451171
poly_degree_fitted_control = 3
poly_degree_polyfit = 4
raw_smooth = 'au_smooth_ma' # 'raw_au'  'au_smooth_ma'
trim_time = 20
signal_chn_name = '_465A'
control_chn_name = '_405A'

# # tmp for mickey data
# signal_chn_name = 'redB'
# control_chn_name = 'isoA'

fpath = os.path.join(r'.\examples\tdt\extracted', 'brandy_03112024_streams_data.feather')
df_data = pd.read_feather(fpath)
streams_path = os.path.join(r'.\examples\tdt\extracted', 'brandy_03112024_streams_info.feather')
streams_info = pd.read_feather(streams_path)
        
# smooth/rolling avg data
df_data = fp_moving_average(df_data, streams_info, 0.1)

# trim time, then separate to channels
df_data['channel_time'] = df_data.groupby('channel').cumcount()/fs
df_data = df_data[df_data['channel_time'] > trim_time].reset_index()

control = df_data[df_data['channel']==control_chn_name][raw_smooth].values
signal = df_data[df_data['channel']==signal_chn_name][raw_smooth].values
time = df_data[df_data['channel']==signal_chn_name]['channel_time'].values

df_processed = pd.DataFrame()

# Generate fitted 405 signal
fit_signal_control = np.polyfit(control, signal, poly_degree_fitted_control)
control_fitted = np.polyval(fit_signal_control, control)

# Generate polynomial fit for signal
fit_poly_signal = np.polyfit(time, signal, poly_degree_polyfit)
poly_signal = np.polyval(fit_poly_signal, time)

# Generate polynomial fit for control
fit_poly_control = np.polyfit(time, control, poly_degree_polyfit)
poly_control = np.polyval(fit_poly_control, time)

# Combine generated data
df_processed['time'] = time
df_processed['control'] = control
df_processed['signal'] = signal

df_processed['control_fitted'] = control_fitted
df_processed['poly_signal'] = poly_signal
df_processed['poly_control'] = poly_control

df_processed['delta_signal_fitted_control'] = signal - control_fitted 
df_processed['delta_signal_poly'] = signal - poly_signal # headfixed
df_processed['delta_control_poly'] = control - poly_control

# Generate z-scores
df_processed['delta_signal_fitted_control_zscore'] = (df_processed['delta_signal_fitted_control'] - df_processed['delta_signal_fitted_control'].mean()) / df_processed['delta_signal_fitted_control'].std() # freely-moving
df_processed['delta_signal_poly_zscore'] = (df_processed['delta_signal_poly'] - df_processed['delta_signal_poly'].mean()) / df_processed['delta_signal_poly'].std()
df_processed['delta_control_poly_zscore'] = (df_processed['delta_control_poly'] - df_processed['delta_control_poly'].mean()) / df_processed['delta_control_poly'].std()

# DFF calculations
mean_100_plus_delta_signal_fitted_control = df_processed['delta_signal_fitted_control'].mean() + 100
df_processed['delta_signal_fitted_control_dff'] = (100 + df_processed['delta_signal_fitted_control'] - mean_100_plus_delta_signal_fitted_control) / mean_100_plus_delta_signal_fitted_control

mean_100_plus_delta_signal_poly = df_processed['delta_signal_poly'].mean() + 100
df_processed['delta_signal_poly_dff'] = (100 + df_processed['delta_signal_poly'] - mean_100_plus_delta_signal_poly) / mean_100_plus_delta_signal_poly

mean_100_plus_delta_control_poly = df_processed['delta_control_poly'].mean() + 100
df_processed['delta_control_poly_dff'] = (100 + df_processed['delta_control_poly'] - mean_100_plus_delta_control_poly) / mean_100_plus_delta_control_poly


df_processed.to_csv('simple_data.csv')

plt.figure()

plt.legend(['delta_signal_fitted_control'])
plt.show()

plt.figure()
plt.plot(df_processed['delta_signal_fitted_control_zscore'])
plt.legend(['delta_signal_fitted_control_zscore'])
plt.show()

######

df_top20 = df_data.groupby('channel', group_keys=False).head(20).reset_index(drop=True)

df_top20.to_csv('2022_06_15_abb07_streams_data.csv')

###############

df_simple = pd.read_csv('simple_data.csv')
fpath_tdt_analysis = pd.read_csv('tdt_analysis_data.csv')
fpath_tdt_analysis = fpath_tdt_analysis.dropna(subset=['time'])
fpath_tdt_analysis = fpath_tdt_analysis[fpath_tdt_analysis['fiber_id'] == 1].reset_index()

plt.figure()
plt.title('Isosbestic')
plt.plot(df_simple['time'], control, alpha=0.8)
plt.plot(fpath_tdt_analysis['time'], fpath_tdt_analysis['control'], alpha=0.5)
plt.legend(['simple', 'tdt_anal'])
plt.xlabel('Time (s)')
plt.ylabel('Raw F')
plt.show()


plt.figure()
plt.title('GCaMP 470')
plt.plot(df_simple['time'], signal, alpha=0.5)
plt.plot(fpath_tdt_analysis['time'], fpath_tdt_analysis['signal'], alpha=0.5)
plt.legend(['Wilson_pipeline', 'Validation'])
plt.xlabel('Time (s)')
plt.ylabel('Raw F')
plt.show()


plt.figure()
plt.title('control fitted')
plt.plot(df_simple['control_fitted'], alpha=0.5)
plt.plot(fpath_tdt_analysis['control_fitted'], alpha=0.5)
plt.legend(['simple', 'tdt_anal'])
plt.show()


plt.figure()
plt.title('delta_signal_fitted_control')
plt.plot(df_simple['time'], df_simple['delta_signal_fitted_control'], alpha=0.5)
plt.plot(fpath_tdt_analysis['time'], fpath_tdt_analysis['delta_signal_fitted_control'], alpha=0.5)
plt.legend(['validation', 'wilson'])
plt.xlabel('Time (s)')
plt.ylabel('Raw F')
plt.show()


########

# testing different models

import numpy as np
import os
import pandas as pd

tmp_data = pd.read_csv('simple_data.csv')

control = tmp_data['control']
signal = tmp_data['signal']

# Generate fitted 405 signal
fit_signal_control = np.polyfit(control, signal, 1)
control_fitted = np.polyval(fit_signal_control, control)

plt.figure()
plt.title('Bad data; LLS-fit')
plt.plot(signal, label='signal')
plt.plot(control, label='control')
plt.plot(control_fitted, label='control_LLS_fitted')
plt.legend()
plt.show()

    
brandy_path = os.path.join(r'C:\Users\stuberadmin\Documents\GitHub\tidy_lab_tools_python\examples\tdt\extracted', 'esr1cre_F04_odor_alb_03112024_streams_session.feather')
brandy_data = pd.read_feather(brandy_path)
control_brandy = brandy_data['control']
signal_brandy = brandy_data['signal']
    
plt.figure()
plt.title('ctrl vs signal')
plt.plot(control-np.mean(control), signal-np.mean(signal), label='bad data')
plt.plot(control_brandy-np.mean(control_brandy), signal_brandy-np.mean(signal_brandy), label='brandy good data', alpha=0.5)
plt.xlabel('control'); plt.ylabel('signal')
plt.legend()
plt.show()

brandy_path = os.path.join(r'C:\Users\stuberadmin\Downloads', 'esr1cre_F04_odor_alb_03112024_streams_data.feather')
brandy_data = pd.read_feather(brandy_path)