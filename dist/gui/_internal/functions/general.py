import pandas as pd
import numpy as np
from typing import List, Optional, Union
import os
import matplotlib.pyplot as plt 

def to_frame(df: Union[pd.DataFrame, pd.Series]) -> pd.DataFrame:
    """
    Convert a pd.Series to a pd.DataFrame.

    Args:
        df (Union[pd.DataFrame, pd.Series]): The input DataFrame or Series.

    Returns:
        pd.DataFrame: The input DataFrame or Series as a DataFrame.
    """
    if isinstance(df, pd.Series):
        df = df.to_frame()
    return df

def plot_dataframe_over_time(df: pd.DataFrame, title=None) -> None:
    """
    Graph a DataFrame over time.

    Args:
        df (pd.DataFrame): The input DataFrame.
        title (str, optional): The title of the graph. Defaults to None.
    """
    for column in df.columns:
        if column != 'time':
            plt.plot(df['time'], df[column], label=column)
    plt.legend()
    plt.title(title)
    plt.xlabel('Time')
    plt.show()

def plot_dataframe(df: pd.DataFrame, title=None, xlabel=None, ylabel=None, x_column=None) -> None:
    """
    Graph a DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame.
        title (str, optional): The title of the graph. Defaults to None.
        xlabel (str, optional): The x-axis label. Defaults to None.
        ylabel (str, optional): The y-axis label. Defaults to None.
        x_column (str, optional): The column to use for the x-axis. Defaults to None.
    """
    df.plot(x=x_column)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

def dir_diff(dir_input: str, dir_output: str, list_suffix: List[str], 
             print_message: int = 0, ignore_suffixes: Optional[List[str]] = None) -> List[str]:
    """
    Return the files located in dir_input that are not already located in dir_output.

    Args:
        dir_input (str): The input directory path.
        dir_output (str): The output directory path.
        list_suffix (List[str]): List of suffixes to remove from file names.
        print_message (int, optional): Whether to print warning messages. Defaults to 0.
        ignore_suffixes (List[str], optional): Suffixes to ignore when comparing. Defaults to None.

    Returns:
        List[str]: List of file names in dir_input but not in dir_output.
    """
    list_dir_input = os.listdir(dir_input)
    list_dir_output = os.listdir(dir_output)

    # Filter out files in dir_output that contain suffixes to ignore
    if ignore_suffixes:
        for ignore_suffix in ignore_suffixes:
            list_dir_output = [fn for fn in list_dir_output if ignore_suffix not in fn]

    for suffix in list_suffix:
        list_dir_input = [fn.replace(suffix, '') for fn in list_dir_input]
        list_dir_output = [fn.replace(suffix, '') for fn in list_dir_output]

    list_dir_input = list(set(list_dir_input))
    list_dir_output = list(set(list_dir_output))

    # If there are any files in the output directory not in the input directory, give a warning
    missing_files = set(list_dir_output) - set(list_dir_input)
    if missing_files and print_message:
        print("WARNING: THE FOLLOWING FILES ARE LOCATED IN dir_output THAT ARE NOT FOUND IN dir_input")
        for fn_missing in missing_files:
            print(fn_missing)
        print("")

    # Create a list of files contained in dir_input but not dir_output
    fns = list(set(list_dir_input) - set(list_dir_output))

    if print_message:
        print(f"Number of files not found in dir_output: {len(fns)}\n")

    return fns

def get_event_bouts(events: pd.DataFrame, event_id_char_of_interest: list, filt_tm: list, filt_n: list, filt_dir: list, id_char: str) -> pd.DataFrame:
    """
    Extract either onset or offsets of bouts of a chosen event based on the number of events prior/post (filt_n) within
    a timeframe prior/post (filt_dir).

    Parameters:
        events (pd.DataFrame): DataFrame from a single session that includes event_id_char (the id of the event) and event_ts (the time of the event).
        event_id_char_of_interest (list): List of event types to filter.
        filt_tm (list): Time for window prior and post event that will be used.
        filt_n (list): Number of events in the window prior and post that will be used.
        filt_dir (list): Logical used to apply to time window ('>', '<', or '==').
        id_char (str): String to replace event_id_char in the output DataFrame.

    Returns:
        pd.DataFrame: DataFrame containing the filtered events.
    """

    # Filter data that match the event of interest
    events_bout_onset: pd.DataFrame = (
        events[events['event_id_char'].isin(event_id_char_of_interest)]
        .assign(n_event_prior=0, n_event_post=0)
        .sort_values(by='event_ts')
    )

    event_tss = events_bout_onset['event_ts'].values

    # For each time stamp, retrieve the number of events preceding and following
    for event in range(events_bout_onset.shape[0]):
        events_bout_onset.at[event, 'n_event_prior'] = np.sum(
            (event_tss > events_bout_onset.at[event, 'event_ts'] - filt_tm[0]) &
            (event_tss < events_bout_onset.at[event, 'event_ts'])
        )

        events_bout_onset.at[event, 'n_event_post'] = np.sum(
            (event_tss > events_bout_onset.at[event, 'event_ts']) &
            (event_tss < events_bout_onset.at[event, 'event_ts'] + filt_tm[1])
        )

    # Filter based on events prior
    if filt_dir[0] == '<':
        events_bout_onset = to_frame(events_bout_onset[events_bout_onset['n_event_prior'] < filt_n[0]])
    elif filt_dir[0] == '>':
        events_bout_onset = to_frame(events_bout_onset[events_bout_onset['n_event_prior'] > filt_n[0]])
    elif filt_dir[0] == '==':
        events_bout_onset = to_frame(events_bout_onset[events_bout_onset['n_event_prior'] == filt_n[0]])

    # Filter based on events post
    if filt_dir[1] == '<':
        events_bout_onset = to_frame(events_bout_onset[events_bout_onset['n_event_post'] < filt_n[1]])
    elif filt_dir[1] == '>':
        events_bout_onset = to_frame(events_bout_onset[events_bout_onset['n_event_post'] > filt_n[1]])
    elif filt_dir[1] == '==':
        events_bout_onset = to_frame(events_bout_onset[events_bout_onset['n_event_post'] == filt_n[1]])

    # Rename event_id_char so it can be combined with the input dataset
    events_bout_onset = (
        events_bout_onset
        .assign(event_id_char=id_char)
        .drop(columns=['n_event_prior', 'n_event_post'])
    )

    return events_bout_onset
