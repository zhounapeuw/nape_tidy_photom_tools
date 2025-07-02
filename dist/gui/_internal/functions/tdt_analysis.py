import pandas as pd
import numpy as np
from functions.general import get_event_bouts, dir_diff, plot_dataframe_over_time, to_frame
import functions.general
import os
import matplotlib.pyplot as plt
from typing import List, Union
from sklearn.linear_model import HuberRegressor
from sklearn.preprocessing import PolynomialFeatures


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


def fp_moving_zscore(streams: pd.DataFrame, pre_window: int) -> pd.DataFrame:
    # returns moving zsocre of delta_signal_poly_zscore in streams
    #
    # inputs
    # - streams (must contain time, fiber_id, signal_wavelength, and delta_signal_poly_zscore)
    # - pre_window number of samples prior for canculating zscore

    temp_moving_zscore = streams[streams['blockname'] != streams['blockname']].loc[:, ['blockname', 'channel_id', 'time']]
    streams_to_process = streams.loc[:, ['blockname', 'fiber_id', 'signal_wavelength']].drop_duplicates()

    for n_stream in range(len(streams_to_process)):
        # Filter data for a particular stream
        temp_stream = streams_to_process.iloc[n_stream].merge(streams, on=['blockname', 'fiber_id', 'signal_wavelength'])
        
        # Calculate z-score for each data point in the stream
        for n_sample in range(pre_window, len(temp_stream)):
            pre_window_data = temp_stream.iloc[(n_sample - pre_window):n_sample]['delta_signal_poly_zscore']
            pre_window_mu = pre_window_data.mean()
            pre_window_sd = pre_window_data.std()

            temp_loop = pd.DataFrame({
                'pre_window_mu': [pre_window_mu],
                'pre_window_sd': [pre_window_sd],
                'blockname': [temp_stream.iloc[n_sample]['blockname']],
                'channel_id': [temp_stream.iloc[n_sample]['channel_id']],
                'time': [temp_stream.iloc[n_sample]['time']]
            })

            temp_moving_zscore = pd.concat([temp_moving_zscore, temp_loop])

    streams = streams.merge(temp_moving_zscore, on=['blockname', 'time', 'channel_id'])
    streams['delta_signal_poly_zscore'] = (streams['delta_signal_poly_zscore'] - streams['pre_window_mu']) / streams['pre_window_sd']

    return streams

def fp_streams_fitted(streams_data: pd.DataFrame, streams_info: pd.DataFrame, signal_id: int, control_id: int, poly_degree_fitted_control: int, poly_degree_polyfit: int, trim_time: float, fit_model: str) -> pd.DataFrame:
    # returns fits and fitted subtractions of data within streams_data.
    #
    # arguments:
    #  streams_data (tibble)- dataframe from *_streams_data.feather exported from python extractor
    #    - fp data stream variable name needs to be raw_au or au_smooth_ma
    #  streams_info dataframe from *_streams_info.feather exported from python extractor
    #  signal_id  (integer)- numerical id for signal  for blockname
    #  control_id (integer)- numerical id for control for blockname
    #  poly_degree_fitted_control (integer)- number of polynomial degrees used for the polynomial fit between control vs signal
    #  poly_degree_polyfit        (integer)- number of polynomial degrees used for the polynomial fit between signal vs time & control vs time
    #  trim_time (seconds) - amount of time to be trimed off start of session
    #
    # output tidy data has observation of time for each fiber within steams_data.
    #  variables in output dataframe
    #      - blockname
    #      - fiber_id       (character, extracted from name of stream within streams_info)
    #      - time           (seconds, calculated from sample number in streams data and fs within streams_info)
    #      raw signals
    #      - signal         (au raw, signal  is defined by signal_id  found within log_fp)
    #      - control        (au raw, control is defined by control_id found within log_fp)
    #      fits
    #      - control_fitted (au, predicted values from polynomial fit of signal vs control, polynomial degree defined by poly_degree_fitted_control)
    #      - poly_signal:   (au, predicted values from polynomial fit of signal  vs time,   polynomial degree defined by poly_degree_polyfit)
    #      - poly_control:  (au, predicted values from polynomial fit of control vs time,   polynomial degree defined by poly_degree_polyfit)
    #      raw fit deltas
    #      - delta_signal_fitted_control (au, signal -  control_fitted)
    #      - delta_signal_poly           (au, signal -  poly_signal)
    #      - delta_control_poly          (au, control - poly_control)

    if sum(streams_data.columns.str.contains('au_smooth_ma')) > 0:
        fp_value = 'au_smooth_ma'
    else:
        fp_value = 'raw_au'


    streams_data = (streams_data[['blockname', 'channel', fp_value]]
            .merge(streams_info[['blockname', 'fs']].drop_duplicates(), on='blockname')
            .groupby(['blockname', 'channel'], group_keys=True, as_index=False)
            .apply(lambda x: x.reset_index(drop=True).assign(time=(np.arange(len(x)) + 1) / x['fs'].iloc[0]))
            .reset_index(drop=True))


    streams_data = fp_identify_fibers(streams_data, 'channel')

    streams_info_fiber_ids = (streams_info[['blockname', 'name']].drop_duplicates()
        .assign(photodetector_id=lambda x: x['name'].str.replace('[0-9]+', '', regex=True))
        .assign(photodetector_id=lambda x: x['photodetector_id'].str.replace('_', '', regex=True))
        .pipe(fp_identify_fibers, 'photodetector_id'))

    fiber_ids = streams_info_fiber_ids['fiber_id'].unique()
    streams_data_processed = None

    for fiber_id_loop in fiber_ids:
        if streams_data['fiber_id'].str.contains(fiber_id_loop).sum() > 0:
            print(f'     > processing fiber number {fiber_id_loop}')

            loop_data: pd.DataFrame = streams_data.loc[(streams_data['fiber_id'].str.contains(fiber_id_loop))]
            loop_data = to_frame(loop_data[['blockname', 'time', 'channel', 'fiber_id', fp_value]])
            loop_data.reset_index(drop=True, inplace=True)

            # temp_end_time = loop_data['time'].reset_index(drop=True).iloc[-1]
            # loop_data['time'] = loop_data['time'] - previous_end_time
            # previous_end_time = temp_end_time

            loop_data = loop_data[loop_data['time'] > trim_time].reset_index()

            # Separate out channel into signal and control columns; reset index for fitting
            time = loop_data.loc[loop_data['channel'].astype(str).str.contains(str(signal_id))]['time'].reset_index(drop=True)
            signal = loop_data.loc[loop_data['channel'].astype(str).str.contains(str(signal_id))][fp_value].reset_index(drop=True)
            control = loop_data.loc[loop_data['channel'].astype(str).str.contains(str(control_id))][fp_value].reset_index(drop=True)

            if fit_model.lower() == 'huber':
                X_poly = PolynomialFeatures(degree=4, include_bias=True).fit_transform(control.values.reshape(-1, 1))
                # Fit robust regression
                model = HuberRegressor().fit(X_poly, signal.values)
                # Predict fitted baseline from control
                control_fitted = model.predict(X_poly)
            else:
                # Generate fitted 405 signal
                fit_signal_control = np.polyfit(control, signal, poly_degree_fitted_control)
                control_fitted = np.polyval(fit_signal_control, control)

            # Generate polynomial fit for signal
            fit_poly_signal = np.polyfit(time, signal, poly_degree_polyfit)
            poly_signal = np.polyval(fit_poly_signal, time)

            # Generate polynomial fit for control
            fit_poly_control = np.polyfit(time, control, poly_degree_polyfit)
            poly_control = np.polyval(fit_poly_control, time)

            tmp_fiber_data = pd.DataFrame({'time': time})
            
            # Combine generated data with loop
            tmp_fiber_data['control'] = control
            tmp_fiber_data['signal'] = signal
            tmp_fiber_data['fiber_id'] = fiber_id_loop
            tmp_fiber_data['control_fitted'] = pd.Series(control_fitted)
            tmp_fiber_data['poly_signal'] = pd.Series(poly_signal)
            tmp_fiber_data['poly_control'] = pd.Series(poly_control)
            tmp_fiber_data['delta_signal_fitted_control'] = pd.Series(signal - control_fitted)
            tmp_fiber_data['delta_signal_poly'] = pd.Series(signal - poly_signal)
            tmp_fiber_data['delta_control_poly'] = pd.Series(control - poly_control)

            # Generate z-scores
            tmp_fiber_data['delta_signal_fitted_control_zscore'] = (tmp_fiber_data['delta_signal_fitted_control'] - tmp_fiber_data['delta_signal_fitted_control'].mean()) / tmp_fiber_data['delta_signal_fitted_control'].std()
            tmp_fiber_data['delta_signal_poly_zscore'] = (tmp_fiber_data['delta_signal_poly'] - tmp_fiber_data['delta_signal_poly'].mean()) / tmp_fiber_data['delta_signal_poly'].std()
            tmp_fiber_data['delta_control_poly_zscore'] = (tmp_fiber_data['delta_control_poly'] - tmp_fiber_data['delta_control_poly'].mean()) / tmp_fiber_data['delta_control_poly'].std()

            # DFF calculations
            mean_100_plus_delta_signal_fitted_control = tmp_fiber_data['delta_signal_fitted_control'].mean() + 100
            tmp_fiber_data['delta_signal_fitted_control_dff'] = (100 + tmp_fiber_data['delta_signal_fitted_control'] - mean_100_plus_delta_signal_fitted_control) / mean_100_plus_delta_signal_fitted_control

            mean_100_plus_delta_signal_poly = tmp_fiber_data['delta_signal_poly'].mean() + 100
            tmp_fiber_data['delta_signal_poly_dff'] = (100 + tmp_fiber_data['delta_signal_poly'] - mean_100_plus_delta_signal_poly) / mean_100_plus_delta_signal_poly

            mean_100_plus_delta_control_poly = tmp_fiber_data['delta_control_poly'].mean() + 100
            tmp_fiber_data['delta_control_poly_dff'] = (100 + tmp_fiber_data['delta_control_poly'] - mean_100_plus_delta_control_poly) / mean_100_plus_delta_control_poly

            tmp_fiber_data['blockname'] = loop_data['blockname'][0] # note this only works if block name is homogenous in loop_data

            if streams_data_processed is None:
                streams_data_processed = tmp_fiber_data
            else:
                streams_data_processed = pd.concat([streams_data_processed, tmp_fiber_data])


    if streams_data_processed is None:
        print("No data to process in fp_streams_fitted function.")
        return pd.DataFrame()

    # return streams_data_processed.reset_index()
    return streams_data_processed
    

def fp_identify_fibers(df: pd.DataFrame, id_name: str) -> pd.DataFrame:
    """
    Return fiber_id for each photodetector present in the dataframe.

    Assumes that photosensor A and photosensor B belong to tdt bank 1 and correspond to fiber 1.
    Assumes that photosensor C and photosensor D belong to tdt bank 1 and correspond to fiber 2.

    :param df: DataFrame containing the photodetector information.
    :param id_name: Name of the column containing photodetector IDs.
    :return: DataFrame with an additional 'fiber_id' column.
    """
    df['fiber_id'] = df[id_name].apply(lambda x: '1' if any(char in x for char in ['A', 'B']) else '2')
    df['fiber_id'] = df['fiber_id'].astype(str)
    return df


def fp_downsample(streams_data_smoothed_fitted: pd.DataFrame, fs: float, down_sampled_freq: float) -> pd.DataFrame:
    # Downsamples all streams / processed streams within streams_data_smoothed_fitted to down_sampled_freq (Hz)
    #
    # Arguments:
    # - streams_data_smoothed_fitted (pd.DataFrame) - produced by fp_streams_fitted
    # - fs (float) - sampling rate (can be returned from streams_info)
    # - down_sampled_freq (float) - resulting frequency of output
    #
    # All columns except blockname and fiber_id will be averaged
    # - Resulting time column will be the maximum time within the sampling period

    tm_bins = np.arange(
        np.round(streams_data_smoothed_fitted['time'].min(), 2),
        np.round(streams_data_smoothed_fitted['time'].max(), 2),
        1 / down_sampled_freq
    )  # Generate time bins and set labels to the end of each window

    streams_data_smoothed_fitted_downsampled = (
        streams_data_smoothed_fitted.assign(
            time=pd.cut(streams_data_smoothed_fitted['time'], bins=tm_bins, labels=tm_bins[1:])
        )
        .groupby(['blockname', 'fiber_id', 'signal_wavelength', 'control_wavelength', 'time'])
        .mean(numeric_only=True)
        .reset_index()
        .assign(time=lambda x: x['time'].astype(float))  # Convert time to numeric
        .groupby(['blockname', 'fiber_id'])
        .filter(lambda x: x['time'].notna().all())
        .query('time < time.max()')  # Remove final bin
        .reset_index(drop=True)
    )

    return streams_data_smoothed_fitted_downsampled

def fp_epocs_to_events(epocs_data: pd.DataFrame, key_file) -> pd.DataFrame:
    # Convert the key_file series into a dataframe with a bunch of columns and 1 row
    key_file = key_file.to_frame().T

    # Labels epochs based on epoch id in log_fp and filters tone
    join_ptc_info = key_file.melt(id_vars=['blockname'], 
                                  value_vars=key_file.filter(like='PtC').columns.to_list() +
                                              key_file.filter(like='PC').columns.to_list(),
                                  var_name='name',
                                  value_name='event_id_char')
    
    join_ptc_info['name'] = join_ptc_info['name'].str.replace('/', '', regex=True).str.replace('Pt', 'P', regex=True)
    join_ptc_info = join_ptc_info.dropna(subset=['event_id_char'])

    epocs_data = to_frame(epocs_data[~epocs_data['name'].str.contains('Cam')])
    epocs_data = to_frame(epocs_data[~epocs_data['name'].str.contains('Tick')])

    epocs_data['name'] = epocs_data['name'].str.replace('/', '', regex=True).str.replace('Pt', 'P', regex=True)
    
    epocs_data = epocs_data.merge(join_ptc_info, on=['blockname', 'name'], how='left')
    epocs_data = epocs_data.rename(columns={'onset': 'event_ts',
                                            'offset': 'event_ts_offset',
                                            'name': 'event_id_ttl'})
    
    epocs_data = to_frame(epocs_data[['blockname', 'event_id_ttl', 'event_id_char', 'event_ts', 'event_ts_offset']])
    
    # Filter out starting impulse across digital inputs upon Arduino start
    epocs_data = filter_first_event(epocs_data)
    
    return epocs_data

def filter_first_event(df: pd.DataFrame) -> pd.DataFrame:
    # Filter out first events that are measured with tdt if there are multiple events detected at arduino startup
    # Only works for sessions with 2 or more TTL inputs

    # Input DataFrame must include the following columns:
    #  - blockname (unique identifier for session)
    #  - event_id_char (event ids)
    #  - event_ts (time of events)

    events_to_filter = df[['blockname', 'event_id_char', 'event_ts']].copy()
    events_to_filter.drop_duplicates(inplace=True)
    events_to_filter['event_ts_round'] = events_to_filter['event_ts'].round(5)
    
    events_to_filter = events_to_filter.groupby('blockname').filter(lambda x: (x['event_ts_round'] == x['event_ts_round'].min()).sum() > 1)
    
    if len(events_to_filter) > 1:
        df = df.merge(events_to_filter.to_frame(), on=['blockname', 'event_id_char', 'event_ts'], how='left', indicator=True)
        df = to_frame(df[df['_merge'] == 'left_only'].drop(columns=['_merge']))
    
    return df

def fp_npm_inputs_to_events(npm_input_data: pd.DataFrame, log_npm: pd.DataFrame) -> pd.DataFrame:

    def f_fp_npm_inputs_to_events(npm_input_data: pd.DataFrame, log_npm: pd.DataFrame) -> pd.DataFrame:
        blockname_process = npm_input_data['blockname_temp'].unique()

        # Labels epochs based on epoch id in log_fp and filters tone
        join_npm_inputs = log_npm.loc[log_npm['blockname'].isin(blockname_process.tolist())]\
            [['input0', 'input1']]\
            .melt(var_name='input_id', value_name='event_id_char', value_vars=['input0', 'input1'])\
            .dropna(subset=['event_id_char'])

        if join_npm_inputs.shape[0] == 0:
            print(f"warning: no npm_inputs defined in log_npm for blockname = {blockname_process[0]}")

        npm_input_data = npm_input_data\
            .merge(join_npm_inputs, left_on='input_id', right_index=True, how='left')\
            [['input_id', 'event_id_char', 'event_ts', 'event_ts_offset']]\
            .assign(event_ts=lambda x: x['event_ts'] / 1000)

        return npm_input_data

    npm_input_data['blockname_temp'] = npm_input_data['blockname']
    grouped = npm_input_data.groupby(['session_name', 'blockname'])
    npm_input_data = grouped.apply(lambda x: f_fp_npm_inputs_to_events(x, log_npm)).reset_index(drop=True).drop_duplicates()

    return npm_input_data

def fp_peri_event_time_histogram(streams_data_smoothed_fitted: pd.DataFrame, 
                                 events_filtered: pd.DataFrame, 
                                 peth_pre: float, 
                                 peth_post: float, 
                                 fs: float, 
                                 down_sampled_freq: float) -> pd.DataFrame:

    streams_data_smoothed_fitted_peth = pd.DataFrame()

    for n_event in range(events_filtered.shape[0]):
        loop_event = events_filtered.iloc[n_event]

        loop_streams = streams_data_smoothed_fitted.loc[
            (streams_data_smoothed_fitted['time'] > loop_event['event_ts'] - peth_pre) &
            (streams_data_smoothed_fitted['time'] < loop_event['event_ts'] + peth_post)
        ]

        loop_streams = loop_streams.merge(loop_event.to_frame().T, on='blockname')
        loop_streams['time_rel'] = loop_streams['time'] - loop_streams['event_ts']

        tm_bins = np.arange(-peth_pre, peth_post, 1 / down_sampled_freq)

        # Resample data
        loop_streams['time_rel'] = pd.cut(loop_streams['time_rel'], tm_bins, labels=tm_bins[1:])

        columns_of_interest = ['session_name', 'blockname', 'fiber_id', 'branch_id', 'signal_wavelength', 'control_wavelength', 'event_id_ttl', 'event_id_char', 'event_ts', 'event_number', 'time_rel']
        existing_columns = [col for col in loop_streams.columns if col in columns_of_interest]
        loop_streams = loop_streams.groupby(existing_columns, as_index=False).mean(numeric_only=True)

        loop_streams['time_rel'] = loop_streams['time_rel'].astype(float)
        loop_streams = loop_streams[loop_streams['time_rel'] < loop_streams['time_rel'].max()]

        if n_event == 0:
            streams_data_smoothed_fitted_peth = loop_streams
        else:
            streams_data_smoothed_fitted_peth = pd.concat([streams_data_smoothed_fitted_peth, loop_streams])

    return streams_data_smoothed_fitted_peth

def resample_peri_event_series(df_series: pd.DataFrame, var_series_ts: str, 
                               time_pre: float, time_post: float, down_sampled_freq: float,
                               vars_grouping: list, vars_to_summarise: list) -> pd.DataFrame:
    # resample data produced from fp_peri_event_time_histogram to set frequency by binning time and summarising
    #
    # df_series (pd.DataFrame): dataframe containing perievent data
    # var_series_ts (str): string of the variable time within df_series that contains the relative time stamps of each sample
    # time_pre (float): time range prior to the event (time rel = 0) that you wish to include in output
    # time_post (float): time range following the event (time rel = 0) that you wish to include in output
    # down_sampled_freq (float, Hz): desired downsampled frequency
    # vars_grouping (list of str): variables you wish to group by (file, subject, trial number, etc.)
    # vars_to_summarise (list of str): variables you wish to compute the mean value of for each time bin
    #
    #  all variables not defined in var_series_ts, vars_grouping, or vars_to_summarise will be dropped.
    
    tm_bins = np.arange(-time_pre, time_post, 1 / down_sampled_freq)  # create time bins

    df_series['rel_time_bin'] = pd.cut(df_series[var_series_ts], tm_bins, labels=False, right=False)  # bin var_series_ts using tm_bins
    df_series['rel_time_bin'] = df_series['rel_time_bin'].apply(lambda x: round(x, 3))  # convert time to numeric
    df_series = to_frame(df_series[df_series['rel_time_bin'] < df_series['rel_time_bin'].max()])  # remove final time bin (often contained in only a subset)

    selected_columns = vars_grouping + ['rel_time_bin'] + vars_to_summarise
    df_series = to_frame(df_series[selected_columns])  # select defined variables

    grouped = df_series.groupby(vars_grouping + ['rel_time_bin'])  # group at all variables except vars_to_summarise
    result = grouped[vars_to_summarise].mean().reset_index()  # summarise to mean of vars_to_summarise within time bins
    return result

def identify_fibers(df: pd.DataFrame, id_name: str) -> pd.DataFrame:
    """
    Return fiber_id for each photodetector present in the dataframe.
    
    Assumes that photosensor A and photosensor B belong to tdt bank 1 and correspond to fiber 1.
    Assumes that photosensor C and photosensor D belong to tdt bank 1 and correspond to fiber 2.
    
    :param df: DataFrame containing the data.
    :param id_name: Name of the column to identify fibers.
    :return: DataFrame with an added 'fiber_id' column.
    """
    df['fiber_id'] = np.where(df[id_name].str.contains('A|B', regex=True), '1', np.nan)
    df['fiber_id'] = np.where(df[id_name].str.contains('C|D', regex=True), '2', df['fiber_id'])
    df['fiber_id'] = df['fiber_id'].astype(str)
    
    return df

def timeseries_resample(df: pd.DataFrame, var_time: str, var_signal: str) -> pd.DataFrame:
    down_sampled_freq = df['down_sampled_freq'].unique().astype(float)
    bin_width = 1000 / down_sampled_freq
    tm_bins = np.arange(0, df[var_time].max() + bin_width, bin_width)

    df_resampled = df.copy()
    df_resampled[var_time] = pd.cut(df_resampled[var_time], list(tm_bins), right=False)
    df_resampled[var_time] = df_resampled[var_time].astype(str).astype(float)
    df_resampled = df_resampled.groupby(df_resampled.columns.difference([var_signal])).agg({var_signal: 'mean'}).reset_index()

    return df_resampled

def timeseries_interpolate_rows(df: pd.DataFrame, var_time: str, var_signal: str) -> pd.DataFrame:
    down_sampled_freq = df['down_sampled_freq'].unique().astype(float)
    bin_width = 1000 / down_sampled_freq
    bin_width_digits = len(str(bin_width).split('.')[1])

    timestamp_seq = np.arange(
        min(df[var_time]),
        max(df[var_time]),
        step=bin_width
    ).round(bin_width_digits)

    interpolated_signal = np.interp(timestamp_seq, df[var_time], df[var_signal])

    filled_data = pd.DataFrame({
        var_time: timestamp_seq,
        var_signal: interpolated_signal
    })
    filled_data[var_time] = filled_data[var_time].round(bin_width_digits)

    return filled_data

def convert_beh_arduino_to_fp_time(df_rec_system_tick: pd.DataFrame, df_beh_events_and_tick: pd.DataFrame) -> pd.DataFrame:
    """
    Convert time for Arduino behavioral events to time for general recording.

    Function converts the time stamps using two clock tickers:
        - behavioral Arduino(s) -> recording system
    Program converts the ticker time of the beh Arduino to the recording system time and then interpolates event times
    based on the time step from the previous ticker (e.g. the recording system time for behavioral Arduino tickers
    recorded by the recording system).

    Importantly, in order to be able to interpolate, the recording system must start recording prior to the onset of the
    behavioral Arduino.

    Args:
        df_rec_system_tick (pd.DataFrame): Dataframe with a time variable named 'ts' and an event_id_char variable
            with the value 'tick_clock_in'.
        df_beh_events_and_tick (pd.DataFrame): Dataframe with a time variable named 'ts' and an event_id_char variable
            with values 'tick_clock_out' and all behavioral events.

    Returns:
        pd.DataFrame: df_beh_transform consisting of df_beh_events_and_tick along with ts_rec_system (derived rec
        system times).
    """

    # clock -> FP
    # prepare fp dataframe for replacing clock ts with fp ts
    ticker_rec = (
        df_rec_system_tick.sort_values('ts')
        .rename(columns={'ts': 'ts_rec_system'})
        .assign(tick_n=lambda x: np.arange(len(x)))
        .assign(event_id_char='tick_clock_out')
    )

    # Replace clock tick ts with FP ts and interpolate beh tick
    ticks = (df_beh_events_and_tick.sort_values('ts') \
                                  .query("event_id_char == 'tick_clock_out'") \
                                  .loc[:, ['ts', 'event_id_char']] \
                                  .rename(columns={'ts': 'ts_tick_previous_beh'}))

    df_beh_transform = (
        df_beh_events_and_tick
        .merge(ticks, on=['event_id_char', 'ts'], how='left')
        .assign(ts_tick_previous_beh=lambda x: x['ts_tick_previous_beh'].fillna(method='ffill'))
        .assign(ts_step=lambda x: x['ts'] - x['ts_tick_previous_beh'])
        .groupby('event_id_char')
        .apply(lambda group: group.assign(tick_n=np.arange(len(group))))
        .reset_index(drop=True)
        .merge(ticker_rec, on=['event_id_char', 'tick_n'], how='left')
        .assign(ts_rec_system=lambda x: x['ts_rec_system'] + x['ts_step'])
        .loc[:, ['event_id_char', 'ts', 'ts_rec_system']]
    )

    return df_beh_transform

def fp_filter_events(events: pd.DataFrame, fp_events_of_interest: list) -> pd.DataFrame:
    """
    Event filtering function.

    Parameters:
        events (pd.DataFrame): DataFrame containing event data.
        fp_events_of_interest (list): List of event types to filter.

    Returns:
        pd.DataFrame: Filtered events.
    """
    # to generate new fileters
    #  1. duplicate a chunk below
    #  2. rename filtered event by changing string in if statement
    #  3. change parameters for get_event_bouts
    #  4. rename filtered event by chaning last string in call

    events_filtered = []

    def in_list(fp_events_of_interest, event):
        # return any(event in x for x in fp_events_of_interest)
        return event in fp_events_of_interest

    if in_list(fp_events_of_interest, 'all'):
        events_filtered = [events]

    if in_list(fp_events_of_interest, 'lick_onset'):
        bout_onset = get_event_bouts(events, ['lick', 'lick01', 'lick02', 'lick03', 'lick04', 'lick05'],
                                     [1, 1], [0, 0], ['==', '>'], 'lick_onset')
        events_filtered.extend(bout_onset)

    if in_list(fp_events_of_interest, 'lick_offset'):
        bout_onset = get_event_bouts(events, ['lick', 'lick01', 'lick02', 'lick03', 'lick04', 'lick05'],
                                     [1, 1], [0, 0], ['>', '=='], 'lick_offset')
        events_filtered.extend(bout_onset)

    if in_list(fp_events_of_interest, 'active_rotation_onset'):
        bout_onset = get_event_bouts(events, ['active_rotation'], [3, 3], [0, 8], ['==', '>'], 'active_rotation_onset')
        events_filtered.extend(bout_onset)

    if in_list(fp_events_of_interest, 'inactive_rotation_onset'):
        bout_onset = get_event_bouts(events, ['inactive_rotation'], [3, 3], [0, 8], ['==', '>'], 'inactive_rotation_onset')
        events_filtered.extend(bout_onset)

    if in_list(fp_events_of_interest, 'active_criteria'):
        bout_onset = events[events['event_id_char'] == 'active_criteria']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'inactive_criteria'):
        bout_onset = events[events['event_id_char'] == 'inactive_criteria']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'active_rotation_criteria'):
        bout_onset = events[events['event_id_char'] == 'active_rotation_criteria']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'inactive_rotation_criteria'):
        bout_onset = events[events['event_id_char'] == 'inactive_rotation_criteria']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'active_rotation_criteria_increment'):
        bout_onset = events[events['event_id_char'] == 'active_rotation_criteria_increment']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'inactive_rotation_criteria_increment'):
        bout_onset = events[events['event_id_char'] == 'inactive_rotation_criteria_increment']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'sol_onset'):
        bout_onset = get_event_bouts(events, ['sol'], [1, 1], [0, 0], ['==', '>'], 'sol_onset')
        events_filtered.extend(bout_onset)

    if in_list(fp_events_of_interest, 'pump_onset'):
        bout_onset = events[events['event_id_char'] == 'sol_onset'].assign(event_id_char='pump_onset')
        events_filtered.extend([bout_onset])

    if in_list(fp_events_of_interest, 'airpuff_onset'):
        bout_onset = events[events['event_id_char'] == 'airpuff_onset']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'cue_led'):
        bout_onset = events[events['event_id_char'] == 'cue_led']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'led_onset'):
        bout_onset = events[events['event_id_char'] == 'led_onset']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'access_period'):
        bout_onset = events[events['event_id_char'] == 'access_period']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'spout_extended'):
        bout_onset = events[events['event_id_char'] == 'spout_extended']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'tone_onset'):
        bout_onset = events[events['event_id_char'].isin(['tone_onset', 'tone'])]
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'shock_onset'):
        bout_onset = events[events['event_id_char'] == 'shock_onset']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'opto_stim'):
        bout_onset = events[events['event_id_char'] == 'opto_stim']
        events_filtered.extend([to_frame(bout_onset)])

    if in_list(fp_events_of_interest, 'bout'):
        bout_onset = (
            events[events['event_id_char'].str.contains('bout')]
            .assign(event_id_char=lambda x: x['event_id_char'] + '_onset')
        )
        bout_offset = (
            events[events['event_id_char'].str.contains('bout')]
            .drop(columns=['event_ts'])
            .assign(event_ts='event_ts_offset', event_id_char=lambda x: x['event_id_char'] + '_offset')
        )
        events_filtered.extend([bout_onset, bout_offset])

    if events_filtered:
        events_filtered = pd.concat(events_filtered)
    else:
        events_filtered = events[events['blockname'] != events['blockname']]

    events_filtered_frame = to_frame(events_filtered)
    return events_filtered_frame

def save_graph(df, dir: str, title=None):
    """
    Saves a DataFrame over time.

    Args:
        df (pd.DataFrame): The input DataFrame.
        title (str, optional): The title of the graph. Defaults to None.
    """

    for column in df.columns:
        if column != 'time':
            plt.plot(df['time'], df[column], label=column)

    plt.gcf().set_size_inches(14, 6)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title(title)
    plt.xlabel('Time')

    plt.tight_layout()
    plt.savefig(f'{dir}/{title}.png')
    plt.clf()

def save_graphs(df: list, dir: str, title: str, titles: list):
    """
    Saves multiple plots onto one image (2x2 grid)
    """

    fig, axs = plt.subplots(2, 2)
    fig.suptitle(title)
    plt.gcf().set_size_inches(28, 12)

    for i, df in enumerate(df):
        for column in df.columns:
            if column != 'time':
                axs[i // 2, i % 2].plot(df['time'], df[column], label=column)
                axs[i // 2, i % 2].set_title(titles[i])
                axs[i // 2, i % 2].legend()

    plt.tight_layout()
    plt.savefig(f'{dir}/{title}.png')
    plt.clf()

def save_streams_data_smoothed_graphs(streams_data_smoothed_fitted_downsampled, output_dir, model_name):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    streams_data_smoothed_fitted_downsampled = streams_data_smoothed_fitted_downsampled.drop(['signal_wavelength', 'control_wavelength', 'blockname'], axis=1)

    fiber_id1 = streams_data_smoothed_fitted_downsampled[streams_data_smoothed_fitted_downsampled['fiber_id'] == '1']
    fiber_id2 = streams_data_smoothed_fitted_downsampled[streams_data_smoothed_fitted_downsampled['fiber_id'] == '2']
    fiber_id1 = fiber_id1.drop(['fiber_id'], axis=1)
    fiber_id2 = fiber_id2.drop(['fiber_id'], axis=1)

    # All traces
    # save_graph(fiber_id1, output_dir, "Fiber 1 All Traces")
    # save_graph(fiber_id2, output_dir, "Fiber 2 All Traces")

    # Poly control and signal
    selected_traces_id1 = fiber_id1[['time', 'signal', 'poly_signal', 'control', 'poly_control']]
    selected_traces_id2 = fiber_id2[['time', 'signal', 'poly_signal', 'control', 'poly_control']]
    # save_graph(selected_traces_id1, output_dir, "Fiber 1 Poly Control and Signal")
    # save_graph(selected_traces_id2, output_dir, "Fiber 2 Poly Control and Signal")

    # Raw signal minus polynomial fit z score
    delta_signal_poly_zscore_id1 = fiber_id1[['time', 'delta_signal_poly_zscore']]
    delta_signal_poly_zscore_id2 = fiber_id2[['time', 'delta_signal_poly_zscore']]
    # save_graph(delta_signal_poly_zscore_id1, output_dir, "Fiber 1 Raw Signal Minus Polynomial Fit Z Score")
    # save_graph(delta_signal_poly_zscore_id2, output_dir, "Fiber 2 Raw Signal Minus Polynomial Fit Z Score")

    # Raw signal and 405 LS fit
    control_fitted_id1 = fiber_id1[['time', 'control', 'signal', 'control_fitted', 'delta_signal_fitted_control']]
    control_fitted_id2 = fiber_id2[['time', 'control', 'signal', 'control_fitted', 'delta_signal_fitted_control']]
    # save_graph(control_fitted_id1, output_dir, "Fiber 1 Raw Signal and 405 LS Fit")
    # save_graph(control_fitted_id2, output_dir, "Fiber 2 Raw Signal and 405 LS Fit")

    # Raw signal minus 405 LA fit z score
    delta_signal_fitted_control_zscore_id1 = fiber_id1[['time', 'delta_signal_fitted_control_zscore']]
    delta_signal_fitted_control_zscore_id2 = fiber_id2[['time', 'delta_signal_fitted_control_zscore']]
    # save_graph(delta_signal_fitted_control_zscore_id1, output_dir, "Fiber 1 Raw Signal Minus 405 LS Fit Z Score")
    # save_graph(delta_signal_fitted_control_zscore_id2, output_dir, "Fiber 2 Raw Signal Minus 405 LS Fit Z Score")

    save_graphs([selected_traces_id1, delta_signal_poly_zscore_id1, control_fitted_id1, delta_signal_fitted_control_zscore_id1], output_dir, "Fiber 1", ["Poly Control and Signal", "Raw Signal Minus Signal Self Poly Z Score", "Raw Signal and 405 LS Fit", f"Raw Signal Minus {model_name} Model LS Fit Z Score"])
    save_graphs([selected_traces_id2, delta_signal_poly_zscore_id2, control_fitted_id2, delta_signal_fitted_control_zscore_id2], output_dir, "Fiber 2", ["Poly Control and Signal", "Raw Signal Minus Signal Self Poly Z Score", "Raw Signal and 405 LS Fit", f"Raw Signal Minus {model_name} Model LS Fit Z Score"])


def fp_preprocess(dir_extracted: str, dir_processed: str, log_fp: pd.DataFrame,
                  manual_experiments: Union[None, List] = None, manual_blocknames: Union[None, List] = None,
                  overwrite: int = 0, dir_extracted_arduino: Union[None, str] = None) -> None:

    # Create dir_processed folder if it doesn't exist
    if not os.path.exists(dir_processed):
        os.makedirs(dir_processed)
    
    if dir_extracted_arduino is None:
        dir_extracted_arduino = None

    extracted_files = os.listdir(dir_extracted)

    # File name suffixes in dir_processed and dir_extracted
    suffixes = ['_streams_peth.feather', '_streams_session.feather', '_events_peth.feather',
                '_epocs_data.feather', '_epocs_info.feather', '_streams_data.feather',
                '_streams_info.feather', '_events.feather', '_streams.feather', '_info.feather',
                '_streams_info.csv']

    # Generate key
    key_files = pd.DataFrame({'blockname': extracted_files})

    if key_files.empty:
        print("WARNING: No files in dir_extracted")
        return

    # Remove suffixes
    for suffix in suffixes:
        key_files['blockname'] = key_files['blockname'].str.replace(suffix, '', regex=True)

    key_files: pd.DataFrame = key_files.drop_duplicates()

    # Remove any files from key_files that aren't in log_fp
    data_in_log_fp = log_fp['subject']
    key_files = to_frame(key_files[key_files['blockname'].str.contains('|'.join(data_in_log_fp), case=False, regex=True)])

    # Join info from log_fp
    key_files = key_files.merge(log_fp, on='blockname', how='left')

    # keep files that are flagged to be included
    key_files = key_files[key_files['include'] == 1]

    # Filter if using manual_blocknames
    if manual_blocknames is not None and len(manual_blocknames) > 0:
        temp = key_files[key_files['blockname'].isin(manual_blocknames)]
        key_files = to_frame(temp)

        if key_files.empty:
            print("No data in dir_extracted matching manual_blocknames")
            print("")
            return

    if manual_experiments is not None and len(manual_experiments) > 0:
        temp = key_files[key_files['experiment'].isin(manual_experiments)]
        key_files = to_frame(temp)

        if key_files.empty:
            print("No data in dir_extracted matching manual_experiments")
            print("")
            return

    # Check to determine if data in dir_extracted are already located in dir_processed
    if overwrite == 0:
        fns_diff = dir_diff(dir_extracted, dir_processed, suffixes, print_message=0,
                            ignore_suffixes=['_events_peth.feather', '_events.feather'])

        temp = key_files[key_files['blockname'].isin(fns_diff)]
        key_files = to_frame(temp)

        # Filter if multi_subject processed
        list_dir_output = os.listdir(dir_processed)
        list_dir_output = [filename.replace(suffix, '') for filename in list_dir_output for suffix in suffixes]
        list_dir_output = list(set(list_dir_output))

        temp = key_files[~key_files['blockname_multi_subject'].isin(list_dir_output)]
        key_files = to_frame(temp)

        # Filter block names produced by multi_subject
        if key_files.empty:
            print("No new data to process")
            print("")
            print("To overwrite data, set overwrite = 1")
            return

    # Bulk processing...
    print(f'bulk processing {key_files.shape[0]} files...')

    for n_file in range(key_files.shape[0]):
        key_file = key_files.iloc[n_file]

        blockname_process = key_file['blockname']
        blockname_process_multi_subject = key_file['blockname_multi_subject']

        if pd.isna(key_file['fit_model']):
            key_file['fit_model'] = 'poly'

        print('_______________________________________________________________________________')
        print('')
        print(f'processing block name: {blockname_process}')

        flag_multi_subject = 0

        bank_process = 0
        if not pd.isna(key_file['multi_config']):
            if key_file['multi_config'] == 'multi_subject':
                flag_multi_subject = 1
                bank_id = None

                if not pd.isna(key_file['fiber_id01']):
                    bank_process = 1
                    bank_id = key_file['fiber_id01']

                if not pd.isna(key_file['fiber_id02']):
                    bank_process = 2
                    bank_id = key_file['fiber_id02']

                if not pd.isna(key_file['fiber_id01']) and not pd.isna(key_file['fiber_id02']):
                    print('WARNING MULTI SUBJECT DEFINED IN BOTH fiber_id01 and fiber_id02')
                    break

                if pd.isna(key_file['fiber_id01']) and pd.isna(key_file['fiber_id02']):
                    print('WARNING MULTI SUBJECT NOT DEFINED IN fiber_id01 and fiber_id02')
                    break

                print(f'processing bank number {bank_process} corresponding to {bank_id}')

        signal_ids = to_frame(key_file.filter(like='signal_id')).melt(var_name='signal_id', value_name='channel_wavelength')
        signal_ids = signal_ids.dropna(subset=['channel_wavelength']).sort_values(by='signal_id')
        signal_ids['control_id'] = key_file['control_id']

        # Import extracted data
        
        streams_path = os.path.join(dir_extracted, f'{blockname_process}_streams_info')
        data_path    = os.path.join(dir_extracted, f'{blockname_process}_streams_data')
        streams_info = load_data(streams_path)
        streams_data = load_data(data_path)

        # Filter stream data and info to bank 1 or 2 if recording multi subject
        if flag_multi_subject:
            if bank_process == 1:
                streams_data = to_frame(streams_data[streams_data['channel'].str.contains('[A-B]')])
                streams_info = to_frame(streams_info[streams_info['name'].str.contains('[A-B]')])
            elif bank_process == 2:
                streams_data = to_frame(streams_data[streams_data['channel'].str.contains('[C-D]')])
                streams_info = to_frame(streams_info[streams_info['name'].str.contains('[C-D]')])

        fs = streams_info['fs'].unique()

        # Smooth data
        print(f' - smoothing data using window width {key_file["smoothing_window"]}')
        streams_data_smoothed: pd.DataFrame = fp_moving_average(streams_data, streams_info, key_file['smoothing_window'])
        
        # Polynomial fitting, least squares fitting, subtractions, zscores, and delta f/f
        print(f' - processing {len(signal_ids)} signal found in session')

        streams_data_smoothed_fitted = pd.DataFrame()
        for _, loop_signal_id in signal_ids.iterrows():
            print(f'   ~ processing signal {loop_signal_id["channel_wavelength"]}')

            # filter streams data to a single signal and control for the loop
            streams_data_loop = streams_data_smoothed[
                (streams_data_smoothed['channel'].str.contains(str(loop_signal_id['channel_wavelength'])) |
                streams_data_smoothed['channel'].str.contains(str(loop_signal_id['control_id'])))
            ]

            streams_data_smoothed_fitted_loop = fp_streams_fitted(streams_data_loop, streams_info,
                            loop_signal_id['channel_wavelength'],
                            key_file['control_id'], key_file['poly_degree_fitted_control'],
                            key_file['poly_degree_polyfit'], key_file['trim_time_start'], key_file['fit_model'])

            streams_data_smoothed_fitted_loop['signal_wavelength'] = loop_signal_id['channel_wavelength']
            streams_data_smoothed_fitted_loop['control_wavelength'] = loop_signal_id['control_id']


            if 'streams_data_smoothed_fitted' not in locals():
                streams_data_smoothed_fitted = streams_data_smoothed_fitted_loop
            else:
                streams_data_smoothed_fitted = pd.concat([streams_data_smoothed_fitted, streams_data_smoothed_fitted_loop])

        if not streams_data_smoothed_fitted.empty:

            # down sample stream data
            print(f' - downsampling data to frequency {key_file["down_sampled_freq"]}')
            streams_data_smoothed_fitted_downsampled = fp_downsample(streams_data_smoothed_fitted, fs[0],
                                                                    key_file['down_sampled_freq'])

            # ------Graphing------
            save_streams_data_smoothed_graphs(streams_data_smoothed_fitted_downsampled.copy(), dir_processed + '/' + 'streams_data_smoothed_graphs/' + blockname_process + '/', key_file['fit_model'])
            # ------Graphing------

            flag_events_tdt_epocs = 0
            flag_events_filtered = 0

            epocs_data = None
            epocs_info = None
            epocs_data_file = os.path.join(dir_extracted, f"{blockname_process}_epocs_data.feather")
            if os.path.exists(epocs_data_file):
                flag_events_tdt_epocs = 1
                epocs_data = pd.read_feather(epocs_data_file)
                epocs_info = pd.read_feather(os.path.join(dir_extracted, f"{blockname_process}_epocs_info.feather"))

            print(" - extracting events and determining events of interest for peth")

            if flag_events_tdt_epocs:
                print("   ~ sourcing events from TDT epocs")

                # get events_of_interest
                events_of_interest = key_file['events_of_interest'].split(";")[0].strip()

                print(f" - using events defined in log_fp: {events_of_interest}")
            
                events = fp_epocs_to_events(epocs_data, key_file)
                events = events.groupby('event_id_char', group_keys=False).apply(lambda x: x.assign(event_number=np.arange(1, len(x) + 1))).reset_index(drop=True)

                if events.empty:
                    print("   * NO EVENTS DETECTED in tdt epoch data *")
                    flag_events_tdt_epocs = 0
                else:
                    # return filtered events based on filters defined in fp_events_of_interest
                    if np.sum(~events['event_id_char'].isna()) > 0:
                        # return empty filtered events if events of interest are listed in log_fp

                        # use clock to align arduino data to fp data
                        if "clock" in events['event_id_char'].unique():
                            print("    ~ sourcing events from arduino file synched to TDT clock")

                            blockname_alignment = blockname_process_multi_subject if flag_multi_subject else blockname_process
                            fn_arduino_events = os.path.join(dir_extracted_arduino, f"{blockname_alignment}_event.csv")

                            if not os.path.exists(fn_arduino_events):
                                print("!error: no extracted arduino file exists")
                                print(fn_arduino_events)
                                return

                            df_beh_events_and_tick = pd.read_csv(fn_arduino_events)
                            df_beh_events_and_tick['event_ts'] = df_beh_events_and_tick['event_ts'] / 1000
                            df_beh_events_and_tick = df_beh_events_and_tick.rename(columns={'event_ts': 'ts'})

                            ticks = df_beh_events_and_tick[df_beh_events_and_tick['event_id_char'] == "tick_clock_out"].assign(tick_number=np.arange(1, len(ticks) + 1))

                            df_rec_system_tick = events[events['event_id_char'] == "clock"].tail(len(ticks)).loc[:, ['event_ts']].rename(columns={'event_ts': 'ts'}).assign(event_id_char="tick_clock_in")

                            events = convert_beh_arduino_to_fp_time(df_rec_system_tick, df_beh_events_and_tick)
                            events = events[~events['ts_rec_system'].isna()].assign(blockname=blockname_process).rename(columns={'ts_rec_system': 'event_ts', 'ts': 'event_ts_arduino'}).groupby('event_id_char').apply(lambda x: x.assign(event_number=np.arange(1, len(x) + 1))).reset_index(drop=True)
                        events_filtered = fp_filter_events(events, events_of_interest)
                        flag_events_filtered = 1
                    else:
                        print("   * NO SPECIFICED EVENTS IN log_fp$events_of_interest *")
            else:
                print("   * no TDT epocs file exists")

            # If there are manual time stamps, read and join with filtered events
            events_manual_fn = key_file['events_manual_fn']  # get events_manual_fn

            if not pd.isna(events_manual_fn):
                print(f"   ~ sourcing events from manual file: {events_manual_fn}")

                events_manual = pd.read_csv(events_manual_fn)
                events_manual = events_manual[events_manual['blockname'] == blockname_process].assign(event_id_ttl="manual")

                if events_manual.empty:
                    print("   * NO EVENTS DETECTED in specified manual file *")
                else:
                    events_manual_filtered = fp_filter_events(events_manual, events_of_interest)

                    if flag_events_filtered:
                        events = pd.concat([events, events_manual], ignore_index=True)
                        events_filtered = pd.concat([events_filtered, events_manual_filtered], ignore_index=True)
                    else:
                        events = events_manual
                        events_filtered = events_manual_filtered
                        flag_events_filtered = 1

                    if events_manual_filtered.empty:
                        print("   * manual events file found but no events made it through filter")

            else:
                print("   * no specified manual events file in log_fp")

            if flag_events_filtered:
                if events_filtered.empty:
                    flag_events_filtered = 0
                    print(" * no perievent time histogram: no filtered time stamps")


            # Perform perievent time histogram if there are events detected
            if flag_events_filtered:
                print(f" - computing peth for {len(events_filtered)} events")

                events_filtered = events_filtered.groupby('event_id_char', group_keys=False).apply(lambda x: x.assign(event_number=np.arange(1, len(x) + 1))).reset_index(drop=True)

                streams_data_smoothed_fitted_peth_downsampled = fp_peri_event_time_histogram(
                    streams_data_smoothed_fitted,
                    events_filtered,
                    key_file['peth_pre'],
                    key_file['peth_post'],
                    fs,
                    key_file['down_sampled_freq']
                )

            # Save processed data
            print(f" - saving files to directory {dir_processed}")

            if flag_multi_subject:
                blockname_process = blockname_process_multi_subject

                streams_data_smoothed_fitted_downsampled = streams_data_smoothed_fitted_downsampled.assign(blockname_multi_subject=blockname_process)
                events = events.assign(blockname_multi_subject=blockname_process)
                events_filtered = events_filtered.assign(blockname_multi_subject=blockname_process)
                streams_data_smoothed_fitted_peth_downsampled = streams_data_smoothed_fitted_peth_downsampled.assign(blockname_multi_subject=blockname_process)

            print(f"   ~ {blockname_process}_streams_session.feather")

            streams_data_smoothed_fitted_downsampled.assign(blockname=blockname_process).to_feather(
                os.path.join(dir_processed, f"{blockname_process}_streams_session.feather")
            )
            
            streams_data_smoothed_fitted_downsampled.to_csv(os.path.join(dir_processed, f"{blockname_process}_streams_session.csv"))

            if flag_events_filtered:  # Save events and peth if they exist
                print(f"   ~ {blockname_process}_events.feather")

                events.assign(blockname=blockname_process).to_feather(
                    os.path.join(dir_processed, f"{blockname_process}_events.feather")
                )

                if len(events_filtered) > 0:
                    print(f"   ~ {blockname_process}_events_peth.feather")

                    events_filtered.assign(blockname=blockname_process).to_feather(
                        os.path.join(dir_processed, f"{blockname_process}_events_peth.feather")
                    )

                    print(f"   ~ {blockname_process}_streams_peth.feather")

                    streams_data_smoothed_fitted_peth_downsampled.assign(blockname=blockname_process).reset_index().to_feather(
                        os.path.join(dir_processed, f"{blockname_process}_streams_peth.feather")
                    )

            
            
            print("")
            print("fp_preprocess() complete")

    return streams_data_smoothed_fitted_downsampled 

if __name__ == "__main__":

    dir_extracted = r'.\examples\tdt\extracted'
    dir_processed = r'.\examples\tdt\processed'
    log_fp = r'examples\tdt\log_data_fp_tdt.csv'
    log_fp = pd.read_csv(log_fp)
    fp_preprocess(dir_extracted, dir_processed, log_fp, overwrite=1)

    # import pandas as pd
    # tmp = pd.read_feather(r"C:\Users\stuberadmin\Documents\GitHub\tidy_lab_tools_python\examples\tdt\extracted\2022_06_15_abb07_streams_info.feather")