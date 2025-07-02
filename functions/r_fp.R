# Title     : r_fp.R
# Objective : To house functions for performing tidy analysis of fiber photometry data produced by the TDT system
# Created by: Adam Gordon-Fennell (agg2248@uw.edu)
#             Garret Stuber Lab, University of Washington
# Created on: 2/15/2022

library(tidyverse)
library(dplyr)
library(magrittr)

# Import functions from r_general.R
# source("r_general.R")
source("C:/Users/stuberadmin/Documents/GitHub/tidy_lab_tools_python/functions/r_general.R")

# preprocessing functions ----------------------------------------------------------------------------------------------
fp_moving_average <- function(streams_data, streams_info, smoothing_window) {
  # returns moving average of raw_au as au_smooth_ma
  #
  # parameters----
  #  streams_data: dataframe from *_streams_data.feather exported from python extractor
  #    - fp data stream variable name needs to be raw_au
  #  streams_info: dataframe from *_streams_info.feather exported from python extractor
  #  smoothing_window: duration of window (seconds) for moving average to be applied
  #    - smoothing window of 0 or NA will return streams_data without applying a moving average

  if (smoothing_window > 0) { # if smoothing...
    require(zoo)

    # extract channel info including sampling rate
    channel_ids <- streams_info %>%
      select(blockname, name, fs) %>%
      unique()

    # calculate moving average for each channel independently
    for (channel_id_index in seq(1, dim(channel_ids)[1])) { # for each channel...

      # filter streams_data to a single channel
      streams_data_loop <- streams_data %>%
        filter(
          blockname == channel_ids$blockname[channel_id_index],
          channel == channel_ids$name[channel_id_index]
        )

      # calculate moving average
      streams_data_loop <- streams_data_loop %>%
        mutate(au_smooth_ma = rollmean(
          x = raw_au, # variable to be smoothed
          k = round(channel_ids$fs[channel_id_index] * smoothing_window), # time window converted to sample number
          align = "center", # align window to center
          na.pad = TRUE
        )) # pad to maintain vector length

      # combine channels into a single dataframe
      if (channel_id_index == 1) {
        streams_data_smoothed <- streams_data_loop
      } else {
        streams_data_smoothed <- streams_data_loop %>% bind_rows(streams_data_smoothed)
      }
    }

    return(streams_data_smoothed %>% filter(!is.na(au_smooth_ma))) # return dataframe with smoothed variable
  } else { # if not smoothing...

    return(streams_data) # return original dataframe without smoothed variable
  }
}

fp_moving_average_grouped <- function(streams_data, log_npm) {
  # returns moving average of raw_au as au_smooth_ma
  #
  # parameters----
  #  streams_data: dataframe from *_streams_data.feather exported from extractor
  #    - fp data stream variable name needs to be raw_au
  #  log_npm
  #    - must include session_name, blockname, and smoothing_window, and down_sampled_freq

  require(zoo)

  # duplicate blockname for group_modify
  streams_data <- streams_data %>% mutate(blockname_temp = blockname)


  # calculate moving average
  f_moving_average <- function(df) {
    k_samples <- df %>%
      mutate(k_samples = down_sampled_freq * smoothing_window) %>%
      pull(k_samples) %>%
      unique()

    if (k_samples > 0) {
      streams_data_smoothed <- df %>%
        mutate(au_smooth_ma = rollmean(
          x = raw_au, # variable to be smoothed
          k = k_samples, # time window converted to sample number
          align = "center", # align window to center
          na.pad = TRUE
        )) # pad to maintain vector length

      streams_data_smoothed <- streams_data_smoothed %>%
        filter(!is.na(au_smooth_ma)) %>%
        select(-down_sampled_freq, -smoothing_window, -blockname_temp)

      return(streams_data_smoothed)
    } else {
      return(df)
    }
  }

  streams_data_smoothed <- streams_data %>%
    left_join(log_npm %>% select(blockname, branch_id, down_sampled_freq, smoothing_window), by = c("blockname", "branch_id")) %>%
    group_by(session_name, blockname, branch_id, channel_wavelength) %>%
    group_modify(~ f_moving_average(.))

  return(streams_data_smoothed)
}



fp_moving_zscore <- function(streams, pre_window) {
  # returns moving zsocre of delta_signal_poly_zscore in streams
  #
  # inputs
  # - streams (must contain time, fiber_id, signal_wavelength, and delta_signal_poly_zscore)
  # - pre_window number of samples prior for canculating zscore

  temp_moving_zscore <- streams %>%
    filter(blockname != blockname) %>%
    select(blockname, channel_id, time)


  streams_to_process <- streams %>%
    select(blockname, fiber_id, signal_wavelength) %>%
    unique()

  for (n_stream in 1:nrow(streams_to_process)) {
    # filter to data to particular stream
    temp_stream <- streams_to_process[n_stream, ] %>%
      left_join(streams, by = c("blockname", "fiber_id", "signal_wavelength"))

    # calculate zscore for each data point in the stream
    for (n_sample in 1 + pre_window:nrow(temp_stream)) {
      temp_loop <- temp_stream[(n_sample - pre_window):(n_sample - 1), ] %>%
        summarise(
          pre_window_mu = mean(delta_signal_poly_zscore),
          pre_window_sd = sd(delta_signal_poly_zscore)
        ) %>%
        mutate(
          blockname = temp_stream$blockname[n_sample],
          channel_id = temp_stream$channel_id[n_sample],
          time = temp_stream$time[n_sample]
        )

      temp_moving_zscore <- temp_moving_zscore %>% bind_rows(temp_loop)
    }
  }

  streams <- streams %>%
    left_join(temp_moving_zscore, by = c("blockname", "time", "channel_id")) %>%
    mutate(delta_signal_poly_zscore_moving = (delta_signal_poly_zscore - pre_window_mu) / pre_window_sd)

  return(streams)
}


fp_streams_fitted <- function(streams_data, streams_info, signal_id, control_id, poly_degree_fitted_control, poly_degree_polyfit, trim_time) {
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
  #
  #


  if (sum(str_detect(colnames(streams_data), "au_smooth_ma")) > 0) {
    fp_value <- "au_smooth_ma"
  } else {
    fp_value <- "raw_au"
  }

  streams_data <- streams_data %>%
    select(blockname, channel, all_of(fp_value)) %>%
    left_join(streams_info %>% select(blockname, fs) %>% unique(), by = "blockname") %>%
    group_by(blockname, channel) %>%
    mutate(time = row_number() / fs) %>%
    ungroup() %>%
    fp_identify_fibers("channel")

  streams_info_fiber_ids <- streams_info %>%
    select(blockname, name) %>%
    unique() %>%
    mutate(photodetector_id = gsub("[[:digit:]]+", "", name)) %>%
    mutate(photodetector_id = gsub("_", "", photodetector_id)) %>%
    fp_identify_fibers("photodetector_id")

  fiber_ids <- unique(streams_info_fiber_ids$fiber_id)

  fiber_loop <- 1

  for (fiber_id_loop in fiber_ids) { # for each fiber id within the blockname
    if (sum(str_detect(streams_data$fiber_id, fiber_id_loop)) > 0) {
      print(str_c("     > processing fiber number ", fiber_id_loop))

      # filter dataset to a single fiber ID and reformat such that each row is a single time point
      loop_data <- streams_data %>%
        filter(str_detect(fiber_id, fiber_id_loop)) %>%
        filter(time > trim_time) %>%
        select(blockname, time, channel, fiber_id, all_of(fp_value)) %>%
        mutate(channel_rename = ifelse(grepl(signal_id, channel), "signal", NA)) %>%
        mutate(channel_rename = ifelse(grepl(control_id, channel), "control", channel_rename)) %>%
        select(-channel) %>%
        spread(channel_rename, fp_value)

      print(loop_data)

      # extract signal and isobesitic vectors
      signal <- loop_data %>% pull(signal)
      control <- loop_data %>% pull(control)
      time <- loop_data %>% pull(time)

      # generate fitted 405 signal
      fit_signal_control <- lm(signal ~ poly(control, poly_degree_fitted_control, raw = TRUE))
      control_fitted <- predict(fit_signal_control, data.frame(x = control))
      
      # generate polynomial fit for signal
      fit_poly_signal <- lm(signal ~ poly(time, poly_degree_polyfit))
      poly_signal <- predict(fit_poly_signal, data.frame(time))

      # generate polynomial fit for control
      fit_poly_control <- lm(control ~ poly(time, poly_degree_polyfit))
      poly_control <- predict(fit_poly_control, data.frame(time))


      # combine generated data with loop
      loop_data <- loop_data %>%
        mutate(
          signal = signal,
          control = control,
          fiber_id = fiber_id_loop,
          control_fitted = control_fitted,
          poly_signal = poly_signal,
          poly_control = poly_control,
          delta_signal_fitted_control = signal - control_fitted,
          delta_signal_poly = signal - poly_signal,
          delta_control_poly = control - poly_control
        )

      # generate zscores
      loop_data <- loop_data %>%
        mutate(
          delta_signal_fitted_control_zscore = (delta_signal_fitted_control - mean(delta_signal_fitted_control)) / sd(delta_signal_fitted_control),
          delta_signal_poly_zscore           = (delta_signal_poly - mean(delta_signal_poly)) / sd(delta_signal_poly),
          delta_control_poly_zscore          = (delta_control_poly - mean(delta_control_poly)) / sd(delta_control_poly),
          delta_signal_fitted_control_dff    = (100 + delta_signal_fitted_control - mean(100 + delta_signal_fitted_control)) / mean(100 + delta_signal_fitted_control),
          delta_signal_poly_dff              = (100 + delta_signal_poly - mean(100 + delta_signal_poly)) / mean(100 + delta_signal_poly),
          delta_control_poly_dff             = (100 + delta_control_poly - mean(100 + delta_control_poly)) / mean(100 + delta_control_poly)
        )

      
      # add in deltas

      if (fiber_loop == 1) {
        streams_data_processed <- loop_data
      } else {
        streams_data_processed <- loop_data %>% bind_rows(streams_data_processed, .)
      }
      
      print(" - printing my_data.csv")
      write.csv(streams_data_processed, "my_data.csv", row.names = FALSE)
      
      fiber_loop <- fiber_loop + 1
    }
  }

  return(streams_data_processed %>% ungroup())
}


fp_normalization_npm <- function(streams_data_smoothed, log_npm) {
  # define function for group_modify
  f_fp_normalization_npm <- function(df, log_npm) {
    # extract signal and isobesitic vectors
    signal <- df %>% pull(signal)
    control <- df %>% pull(control)
    time <- df %>% pull(time)



    blockname_function <- df$blockname_temp %>% unique()

    poly_degree_fitted_control <- log_npm %>%
      filter(blockname == blockname_function) %>%
      pull(poly_degree_fitted_control) %>%
      unique()
    poly_degree_polyfit <- log_npm %>%
      filter(blockname == blockname_function) %>%
      pull(poly_degree_polyfit) %>%
      unique()

    # generate fitted 405 signal
    fit_signal_control <- lm(signal ~ poly(control, poly_degree_fitted_control, raw = TRUE))
    control_fitted <- predict(fit_signal_control, data.frame(x = control))
    
    # generate polynomial fit for signal
    fit_poly_signal <- lm(signal ~ poly(time, poly_degree_polyfit))
    poly_signal <- predict(fit_poly_signal, data.frame(time))

    # generate polynomial fit for control
    fit_poly_control <- lm(control ~ poly(time, poly_degree_polyfit))
    poly_control <- predict(fit_poly_control, data.frame(time))


    # combine generated data with loop
    df <- df %>%
      mutate(
        signal = signal,
        control = control,
        control_fitted = control_fitted,
        poly_signal = poly_signal,
        poly_control = poly_control,
        delta_signal_fitted_control = signal - control_fitted,
        delta_signal_poly = signal - poly_signal,
        delta_control_poly = control - poly_control
      )

    df <- df %>%
      mutate(
        delta_signal_fitted_control_zscore = (delta_signal_fitted_control - mean(delta_signal_fitted_control)) / sd(delta_signal_fitted_control),
        delta_signal_poly_zscore           = (delta_signal_poly - mean(delta_signal_poly)) / sd(delta_signal_poly),
        delta_control_poly_zscore          = (delta_control_poly - mean(delta_control_poly)) / sd(delta_control_poly),
        delta_signal_fitted_control_dff    = (100 + delta_signal_fitted_control - mean(100 + delta_signal_fitted_control)) / mean(100 + delta_signal_fitted_control),
        delta_signal_poly_dff              = (100 + delta_signal_poly - mean(100 + delta_signal_poly)) / mean(100 + delta_signal_poly),
        delta_control_poly_dff             = (100 + delta_control_poly - mean(100 + delta_control_poly)) / mean(100 + delta_control_poly)
      )
    

    
    return(df)
  }


  streams_data_smoothed_fitted <- log_npm %>%
    select(session_name, blockname, branch_id, signal_id01, signal_id02, control_id) %>%
    gather("channel_id", "channel_wavelength", signal_id01:control_id) %>%
    mutate(channel_id = ifelse(channel_id == "control_id", "control", "signal")) %>%
    mutate(channel_wavelength = channel_wavelength %>% as.character()) %>%
    left_join(streams_data_smoothed, ., by = c("session_name", "blockname", "branch_id", "channel_wavelength"))

  streams_data_smoothed_fitted <- streams_data_smoothed_fitted %>%
    arrange(session_name, blockname, branch_id, time) %>%
    select(session_name, blockname, branch_id, time, channel_wavelength, channel_id, au_smooth_ma) %>%
    group_by(session_name, blockname, branch_id, channel_wavelength) %>%
    mutate(cycle_number = row_number()) %>%
    group_by(session_name, blockname, branch_id, cycle_number) %>%
    mutate(time = mean(time)) %>%
    spread(channel_id, au_smooth_ma) %>%
    mutate(control = mean(control, na.rm = "true")) %>%
    filter(!is.na(signal)) %>%
    ungroup()


  streams_data_smoothed_fitted <- streams_data_smoothed_fitted %>%
    group_by(session_name, blockname, branch_id, channel_wavelength) %>%
    mutate(blockname_temp = blockname) %>%
    group_modify(~ f_fp_normalization_npm(., log_npm))

  return(streams_data_smoothed_fitted)
}

fp_downsample <- function(streams_data_smoothed_fitted, fs, down_sampled_freq) {
  # downsamples all streams / processed streams within streams_data_smoothed_fitted to down_sampled_freq (Hz)
  #
  # arguments:
  #  - streams_data_smoothed_fitted (tibble)- produced by fp_streams_fitted
  #  - fs (Hz)- sampling rate (can be returned from streams_info)
  #  - down_sampled_freq (Hz)- resulting frequency of output
  #
  # all columns except blockname and fiber_id will be averaged
  #  - resulting time column will be the maximum time within the sampling period

  tm_bins <- seq(
    min(streams_data_smoothed_fitted$time) %>% round(2),
    max(streams_data_smoothed_fitted$time) %>% round(2),
    1 / down_sampled_freq
  ) # generate time bins and set labels to end of each window

  streams_data_smoothed_fitted_downsampled <- streams_data_smoothed_fitted %>%
    mutate(time = cut(time, tm_bins, tm_bins[2:length(tm_bins)])) %>%
    group_by(blockname, fiber_id, signal_wavelength, control_wavelength, time) %>%
    summarise_all(mean) %>% # summarise all variables using mean
    mutate(time = time %>% as.character() %>% as.double()) %>% # convert time to numeric
    group_by(blockname, fiber_id) %>%
    filter(!is.na(time)) %>%
    filter(time < max(time)) # remove final bin


  return(streams_data_smoothed_fitted_downsampled %>% ungroup())
}

fp_epocs_to_events <- function(epocs_data, key_file) {
  # labels epochs based on epoch id in log_fp and filteres tone
  join_ptc_info <- key_file %>%
    select(blockname, starts_with(c("PtC", "PC"))) %>%
    gather("name", "event_id_char", starts_with(c("PtC", "PC"))) %>%
    mutate(name = name %>% str_remove("/")) %>%
    mutate(name = name %>% str_replace("Pt", "P")) %>%
    filter(!is.na(event_id_char)) # filter to only events defined in event_id_char

  epocs_data <- epocs_data %>%
    filter(!grepl("Cam", name)) %>% # remove cam
    filter(!grepl("Tick", name)) %>% # remove tick from epocs
    mutate(name = name %>% str_remove("/")) %>%
    mutate(name = name %>% str_replace("Pt", "P")) %>%
    left_join(join_ptc_info, by = c("blockname", "name")) %>%
    rename(
      event_ts = onset,
      event_ts_offset = offset,
      event_id_ttl = name
    ) %>%
    select(blockname, event_id_ttl, event_id_char, event_ts, event_ts_offset)

  # filter out starting impulse across digital inputs upon arduino start
  epocs_data <- epocs_data %>%
    filter_first_event()

  # epocs_data <- epocs_data %>%
  # filter(!is.na(event_id_char))

  return(epocs_data)
}

fp_npm_inputs_to_events <- function(npm_input_data, log_npm) {
  f_fp_npm_inputs_to_events <- function(npm_input_data, log_npm) {
    blockname_process <- npm_input_data %>%
      pull(blockname_temp) %>%
      unique()

    # labels epochs based on epoch id in log_fp and filteres tone
    join_npm_inputs <- log_npm %>%
      filter(blockname == blockname_process) %>%
      select(input0, input1) %>%
      gather("input_id", "event_id_char", starts_with("input")) %>%
      filter(!is.na(event_id_char)) # filter to only events defined in event_id_char

    if (nrow(join_npm_inputs) == 0) {
      print(str_c("warning: no npm_inputs defined in log_npm for blockname = ", blockname_process))
    }

    npm_input_data <- npm_input_data %>%
      left_join(join_npm_inputs, by = "input_id") %>%
      select(input_id, event_id_char, event_ts, event_ts_offset) %>%
      mutate(event_ts = event_ts / 1000)

    return(npm_input_data)
  }


  npm_input_data <- npm_input_data %>%
    mutate(blockname_temp = blockname) %>%
    group_by(session_name, blockname) %>%
    group_modify(~ f_fp_npm_inputs_to_events(., log_npm)) %>%
    ungroup() %>%
    unique()

  return(npm_input_data)
}

fp_peri_event_time_histogram <- function(streams_data_smoothed_fitted, events_filtered, peth_pre, peth_post, fs, down_sampled_freq) {
  for (n_event in seq(1, dim(events_filtered)[1])) {
    loop_event <- events_filtered[n_event, ]

    loop_streams <- streams_data_smoothed_fitted %>%
      filter(time > loop_event$event_ts - peth_pre, time < loop_event$event_ts + peth_post) %>%
      left_join(loop_event, by = "blockname") %>%
      mutate(time_rel = time - event_ts)

    tm_bins <- seq(-peth_pre, peth_post, 1 / down_sampled_freq) # generate time bins and set labels to end of each window

    # resample data
    loop_streams <- loop_streams %>%
      mutate(time_rel = cut(time_rel, tm_bins, tm_bins[2:length(tm_bins)])) %>%
      group_by(across(any_of(c("session_name", "blockname", "fiber_id", "branch_id", "signal_wavelength", "control_wavelength", "event_id_ttl", "event_id_char", "event_ts", "event_number", "time_rel")))) %>%
      summarise_all(mean, .groups = "drop") %>% # summarise all variables using mean
      mutate(time_rel = time_rel %>% as.character() %>% as.double()) %>% # convert time to numeric
      filter(time_rel < max(time_rel))

    if (n_event == 1) {
      streams_data_smoothed_fitted_peth <- loop_streams
    } else {
      streams_data_smoothed_fitted_peth <- loop_streams %>% bind_rows(streams_data_smoothed_fitted_peth, .)
    }
  }

  return(streams_data_smoothed_fitted_peth)
}

resample_peri_event_series <- function(df_series, var_series_ts, time_pre, time_post, down_sampled_freq, vars_grouping, vars_to_summarise) {
  # resample data produced from fp_peri_event_time_histogram to set frequency by binning time and summarising
  #
  # df_series (dataframe) dataframe containing perievent data
  # var_series_ts (string): string of the variable time within df_series that contains the relative time stamps of each sample
  # time_pre (double): time range prior to the the event (time rel = 0) that you wish to include in output
  # time_post (double): time range following to the the event (time rel = 0) that you wish to include in output
  # down_sampled_freq (double, Hz): desired downsampled frequency
  # vars_grouping (vector of strings): variables you wish to group by (file, subject, trial number, etc.)
  # vars_to_summarise (vector of strings): variables you wish to compute the mean value of for each time bin
  #
  #  all variables not defined in var_series_ts, vars_grouping, or vars_to_summarise will be dropped.

  tm_bins <- seq(-time_pre, time_post, 1 / down_sampled_freq) # create time bins

  df_series %>%
    mutate(rel_time_bin = cut(!!as.name(var_series_ts), tm_bins, tm_bins[2:length(tm_bins)])) %>% # bin var_series_ts using tm_bins
    mutate(rel_time_bin = rel_time_bin %>% as.character() %>% as.double() %>% round(3)) %>% # convert time to numeric
    filter(rel_time_bin < max(rel_time_bin)) %>% # remove final time bin (often contained in only a subset)
    select(one_of(vars_grouping, "rel_time_bin", vars_to_summarise)) %>% # select to defined variables
    group_by_at(c(vars_grouping, "rel_time_bin")) %>% # group at all variables except vars_to_summarise
    summarise_all(mean, .groups = "drop") # summarise to mean of vars_to_summarise within time bins
}

fp_identify_fibers <- function(df, id_name) {
  # return fiber_id for each photodetector present in the dataframe
  #
  # assuems that photosensor A and photosensor B belong to tdt bank 1 and corresponds to fiber 1
  # assuems that photosensor C and photosensor D belong to tdt bank 1 and corresponds to fiber 2

  df <- df %>%
    mutate(fiber_id = ifelse(!!as.name(id_name) %>% str_detect("A|B"), 1, NA)) %>%
    mutate(fiber_id = ifelse(!!as.name(id_name) %>% str_detect("C|D"), 2, fiber_id)) %>%
    mutate(fiber_id = as.character(fiber_id))

  return(df)
}

# resample dataset
timeseries_resample <- function(df, var_time, var_signal) {
  bin_width <- 1000 / (df %>% pull(down_sampled_freq) %>% unique())
  tm_bins <- seq(0, max(df %>% pull(!!as.name(var_time))) + bin_width, bin_width)

  df_resampled <- df %>%
    mutate(!!as.name(var_time) := cut(!!as.name(var_time), tm_bins, tm_bins[2:length(tm_bins)])) %>%
    mutate(!!as.name(var_time) := !!as.name(var_time) %>%
      as.character() %>%
      as.double()) %>%
    group_by(across(-!!as.name(var_signal))) %>%
    summarise(!!as.name(var_signal) := mean(!!as.name(var_signal)), .groups = "drop")

  return(df_resampled)
}


# interpolate
timeseries_interpolate_rows <- function(df, var_time, var_signal) {
  bin_width <- 1000 / (df %>% pull(down_sampled_freq) %>% unique())
  bin_width_digits <- nchar(sub(".*\\.", "", as.character(bin_width)))

  timestamp_seq <- seq(
    min(df[[var_time]]),
    max(df[[var_time]]),
    by = bin_width
  ) %>%
    round(digits = bin_width_digits)


  filled_data <- tibble(
    !!as.name(var_time) := timestamp_seq,
    !!as.name(var_signal) := approx(df[[var_time]], df[[var_signal]], xout = timestamp_seq)$y
  ) %>%
    mutate(!!as.name(var_time) := !!as.name(var_time) %>% round(bin_width_digits))

  return(filled_data)
}

convert_beh_arduino_to_fp_time <- function(df_rec_system_tick, df_beh_events_and_tick) {
  #' Convert time for arduino behavioral events to time for general recording
  #'
  #' Function converts the time stamps using two clock tickers:
  #'   - behavioral arduino(s) -> recording system
  #' Program converts the ticker time of the beh arduino to the recording system time and then interpolates event times based on the time step from the previous ticker (e.g. the recording system time for behavioral arduino tickers recorded by the rec ording system)
  #'
  #' Importantly, in order to be able to interpolate, the recording system must start recording prior to the onset of the behavioral arduino
  #'
  #'
  #' @param df_rec_system_tick Dataframe with a time variable named 'ts' and an event_id_char variable with the value 'tick_clock_in'.
  #' @param  df_beh_events_and_tick Dataframe with a time variable named 'ts' and an event_id_char variable with values 'tick_clock_out' and all behavioral events.
  #' @return df_beh_transform consisting of df_beh_events_and_tick along with ts_rec_system (derived rec system times)

  # clock -> FP
  ## prepare fp dataframe for replacing clock ts with fp ts
  ticker_rec <- df_rec_system_tick %>%
    arrange(ts) %>%
    rename(ts_rec_system = ts) %>%
    mutate(tick_n = row_number()) %>%
    mutate(event_id_char = "tick_clock_out")

  ## replace clock tick ts with fp ts and interpolate beh tick
  ticks <- df_beh_events_and_tick %>%
    arrange(ts) %>%
    filter(event_id_char == "tick_clock_out") %>%
    select(ts, event_id_char) %>%
    mutate(ts_tick_previous_beh = ts)

  df_beh_transform <- df_beh_events_and_tick %>%
    left_join(ticks, by = c("event_id_char", "ts")) %>%
    fill(ts_tick_previous_beh) %>%
    mutate(ts_step = ts - ts_tick_previous_beh) %>% # save time step to fill in tick_beh_in
    group_by(event_id_char) %>%
    mutate(tick_n = row_number()) %>%
    ungroup() %>%
    left_join(ticker_rec, by = c("event_id_char", "tick_n")) %>%
    fill(ts_rec_system) %>%
    mutate(ts_rec_system = ts_rec_system + ts_step)

  df_beh_transform <- df_beh_transform %>%
    select(event_id_char, ts, ts_rec_system)

  return(df_beh_transform)
}



# event filtering -----------------------------------------------------------------------------------------------------
fp_filter_events <- function(events, fp_events_of_interest) {
  #
  # to generate new fileters
  #  1. duplicate a chunk below
  #  2. rename filtered event by changing string in if statement
  #  3. change parameters for get_event_bouts
  #  4. rename filtered event by chaning last string in call

  events_filtered <- list()

  if (sum(str_detect(fp_events_of_interest, "all")) > 0) {
    events_filtered <- list(events)
  }

  if (sum(str_detect(fp_events_of_interest, "lick_onset")) > 0) {
    bout_onset <- get_event_bouts(events, c("lick", "lick01", "lick02", "lick03", "lick04", "lick05"), c(1, 1), c(0, 0), c("==", ">"), "lick_onset")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "lick_offset")) > 0) {
    bout_onset <- get_event_bouts(events, c("lick", "lick01", "lick02", "lick03", "lick04", "lick05"), c(1, 1), c(0, 0), c(">", "=="), "lick_offset")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "active_rotation_onset")) > 0) {
    bout_onset <- get_event_bouts(events, c("active_rotation"), c(3, 3), c(0, 8), c("==", ">"), "active_rotation_onset")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "inactive_rotation_onset")) > 0) {
    bout_onset <- get_event_bouts(events, c("inactive_rotation"), c(3, 3), c(0, 8), c("==", ">"), "inactive_rotation_onset")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "active_criteria")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "active_criteria")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "inactive_criteria")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "inactive_criteria")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "active_rotation_criteria")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "active_rotation_criteria")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "inactive_rotation_criteria")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "inactive_rotation_criteria")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "active_rotation_criteria_increment")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "active_rotation_criteria_increment")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "inactive_rotation_criteria_increment")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "inactive_rotation_criteria_increment")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "sol_onset")) > 0) {
    bout_onset <- get_event_bouts(events, c("sol"), c(1, 1), c(0, 0), c("==", ">"), "sol_onset")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "pump_onset")) > 0) {
    bout_onset <- events %>%
      filter(event_id_char == "sol_onset") %>%
      mutate(event_id_char = "pump_onset")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "airpuff_onset")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "airpuff_onset")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "cue_led")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "cue_led")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "led_onset")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "led_onset")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "access_period")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "access_period")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "spout_extended")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "spout_extended")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "tone_onset")) > 0) {
    bout_onset <- events %>% filter(event_id_char %in% c("tone_onset", "tone"))
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "shock_onset")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "shock_onset")
    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (sum(str_detect(fp_events_of_interest, "opto_stim")) > 0) {
    bout_onset <- events %>% filter(event_id_char == "opto_stim")
    events_filtered <- c(events_filtered, list(bout_onset))
  }


  if (sum(str_detect(fp_events_of_interest, "bout")) > 0) {
    bout_onset <- events %>%
      filter(str_detect(event_id_char, "bout")) %>%
      mutate(event_id_char = event_id_char %>% str_c("_onset"))

    bout_onset <- events %>%
      filter(str_detect(event_id_char, "bout")) %>%
      select(-event_ts) %>%
      rename(event_ts = event_ts_offset) %>%
      mutate(event_id_char = event_id_char %>% str_c("_offset")) %>%
      bind_rows(bout_onset, .)

    events_filtered <- c(events_filtered, list(bout_onset))
  }

  if (length(events_filtered) >= 1) {
    events_filtered <- bind_rows(events_filtered)
  } else {
    events_filtered <- events %>% filter(blockname != blockname)
  }

  return(events_filtered)
}

filter_first_event <- function(df) {
  # filter out first events that are measured with tdt if there are multiple events detected at arduino startup
  # * only works for sessions with 2 or more TTL inputs
  #
  # input dataframe must include the following variables
  #  - blockname (unique identifier for session)
  #  - event_id_char (event ids)
  #  - event_ts (time of events)

  events_to_filter <- df %>%
    ungroup() %>%
    select(blockname, event_id_char, event_ts) %>%
    unique() %>%
    group_by(blockname) %>%
    mutate(event_ts_round = round(event_ts, 5)) %>%
    filter(event_ts_round == min(event_ts_round))

  if (nrow(events_to_filter) > 1) {
    df <- df %>%
      anti_join(events_to_filter)
  }

  return(df)
}

# pre processing executions---------------------------------------------------------------------------------------------

fp_preprocess <- function(dir_extracted, dir_processed, log_fp, manual_experiments = NA, manual_blocknames = NA, overwrite = 0, dir_extracted_arduino = NA) {
  if (missing(dir_extracted_arduino)) {
    dir_extracted_arduino <- NA
  }

  extracted_files <- list.files(dir_extracted)

  # file name suffixes in dir_processed and dir_extracted
  suffixes <- c(
    "_streams_peth.feather",
    "_streams_session.feather",
    "_events_peth.feather",
    "_epocs_data.feather",
    "_epocs_info.feather",
    "_streams_data.feather",
    "_streams_info.feather",
    "_events.feather",
    "_streams.feather",
    "_info.feather"
  )

  # generate key
  key_files <- tibble("blockname" = extracted_files)

  if (nrow(key_files) == 0) {
    print("WARNING: No files in dir_extracted")
    return()
  }

  # remove suffixes
  for (suffix in suffixes) {
    key_files <- key_files %>% mutate(blockname = blockname %>% str_remove(suffix))
  }

  key_files <- key_files %>% unique()

  # join info from log_fp
  key_files <- key_files %>%
    left_join(log_fp, by = "blockname")

  # filter if using manual_blocknames
  if (sum(!is.na(manual_blocknames)) > 0) {
    key_files <- key_files %>%
      filter(blockname %in% manual_blocknames)

    if (nrow(key_files) == 0) {
      print("No data in dir_extracted matching manual_blocknames")
      print("")

      return()
    }
  }

  if (sum(!is.na(manual_experiments)) > 0) {
    key_files <- key_files %>%
      filter(experiment %in% manual_experiments)

    if (nrow(key_files) == 0) {
      print("No data in dir_extracted matching manual_experiments")
      print("")

      return()
    }
  }

  # check to determine if data in dir_extracted are already located in dir_processed
  if (overwrite == 0) {
    fns_diff <- dir_diff(dir_extracted, dir_processed, suffixes, print_message = 0, ignore_suffixes = c("_events_peth.feather", "_events.feather"))

    key_files <- key_files %>%
      filter(blockname %in% fns_diff)

    # filter if multi_subject processed
    list_dir_output <- list.files(dir_processed)

    for (suffix in suffixes) {
      list_dir_output <- list_dir_output %>%
        str_remove(suffix)
    }

    list_dir_output <- list_dir_output %>% unique()

    key_files <- key_files %>%
      filter(!blockname_multi_subject %in% list_dir_output)

    # filter block names produced by multi_subject

    if (nrow(key_files) == 0) {
      print("No new data to process")
      print("")
      print("to overwrite data, set overwrite = 1")

      return()
    }
  }


  # Bulk processing...
  print(str_c("bulk processing ", dim(key_files)[1], " files..."))

  for (n_file in seq(1, dim(key_files)[1])) {
    key_file <- key_files[n_file, ]

    blockname_process <- key_file$blockname
    blockname_process_multi_subject <- key_file$blockname_multi_subject # save output filename for multi_subject

    print("_______________________________________________________________________________")
    print("")
    print(str_c("processing block name: ", blockname_process))

    flag_multi_subject <- 0

    if (!is.na(key_file$multi_config)) {
      if (key_file$multi_config == "multi_subject") {
        flag_multi_subject <- 1

        if (!is.na(key_file$fiber_id01)) {
          bank_process <- 1
          bank_id <- key_file$fiber_id01
        }

        if (!is.na(key_file$fiber_id02)) {
          bank_process <- 2
          bank_id <- key_file$fiber_id02
        }

        if (!is.na(key_file$fiber_id01) & !is.na(key_file$fiber_id02)) {
          print(str_c("WARNING MULTI SUBJECT DEFINED IN BOTH fiber_id01 and fiber_id02"))
          break
        }

        if (is.na(key_file$fiber_id01) & is.na(key_file$fiber_id02)) {
          print(str_c("WARNING MULTI SUBJECT NOT DEFINED IN fiber_id01 and fiber_id02"))
          break
        }

        print(str_c("processing bank number ", bank_process, " corresponding to ", bank_id))
      }
    }

    # determine signals within blockname
    signal_ids <- key_file %>%
      select(blockname, starts_with("signal_id")) %>%
      gather("signal_id", "channel_wavelength", starts_with("signal_id")) %>%
      filter(!is.na(channel_wavelength)) %>%
      arrange(signal_id) %>%
      mutate(control_id = key_file$control_id)

    # import extracted  data
    streams_info <- arrow::read_feather(str_c(dir_extracted, "/", blockname_process, "_streams_info.feather"))
    streams_data <- arrow::read_feather(str_c(dir_extracted, "/", blockname_process, "_streams_data.feather"))

    # filter stream data and info to bank 1 or 2 if recording multi subject
    if (flag_multi_subject) {
      if (bank_process == 1) {
        streams_data <- streams_data %>%
          filter(str_detect(channel, "A") | str_detect(channel, "B"))

        streams_info <- streams_info %>%
          filter(str_detect(name, "A") | str_detect(name, "B"))
      }

      if (bank_process == 2) {
        streams_data <- streams_data %>%
          filter(str_detect(channel, "C") | str_detect(channel, "D"))

        streams_info <- streams_info %>%
          filter(str_detect(name, "C") | str_detect(name, "D"))
      }
    }

    fs <- streams_info %>%
      pull(fs) %>%
      unique()

    # smooth data
    print(str_c(" - smoothing data using window width ", key_file$smoothing_window))

    streams_data_smoothed <-
      fp_moving_average(
        streams_data,
        streams_info,
        key_file$smoothing_window
      )

    # polynomial fitting, least squares fitting, subtractions, zscores, and delta f/f
    print(str_c(" - processing ", nrow(signal_ids), " signal found in session"))

    for (n_signal in 1:nrow(signal_ids)) {
      loop_signal_id <- signal_ids[n_signal, ]

      print(str_c("   ~ processing signal ", loop_signal_id$channel_wavelength))


      # filter streams data to a single signal and control for the loop
      streams_data_loop <- streams_data_smoothed %>%
        filter(
          channel %>% str_detect(loop_signal_id$channel_wavelength %>% as.character()) |
            channel %>% str_detect(loop_signal_id$control_id %>% as.character())
        )

      # perform fits and normaizations
      streams_data_smoothed_fitted_loop <-
        fp_streams_fitted(
          streams_data = streams_data_loop,
          streams_info = streams_info,
          signal_id = loop_signal_id$channel_wavelength %>% as.character(),
          control_id = key_file$control_id %>% as.character(),
          poly_degree_fitted_control = key_file$poly_degree_fitted_control,
          poly_degree_polyfit = key_file$poly_degree_polyfit,
          trim_time = key_file$trim_time_start
        ) %>%
        mutate(
          signal_wavelength = loop_signal_id$channel_wavelength,
          control_wavelength = loop_signal_id$control_id
        )

      if (n_signal == 1) {
        streams_data_smoothed_fitted <- streams_data_smoothed_fitted_loop
      } else {
        streams_data_smoothed_fitted <- streams_data_smoothed_fitted_loop %>%
          bind_rows(streams_data_smoothed_fitted)
      }
    }

    # down sample stream data
    print(str_c(" - downsampling data to frequency ", key_file$down_sampled_freq))

    streams_data_smoothed_fitted_downsampled <-
      fp_downsample(
        streams_data_smoothed_fitted,
        fs,
        key_file$down_sampled_freq
      )

    # extract events-----------------------------------------------------------------------------------------------
    flag_events_tdt_epocs <- 0
    flag_events_filtered <- 0

    # read in epoc data if it exists
    if (file.exists(str_c(dir_extracted, "/", blockname_process, "_epocs_data.feather"))) {
      flag_events_tdt_epocs <- 1
      epocs_data <- arrow::read_feather(str_c(dir_extracted, "/", blockname_process, "_epocs_data.feather"))
      epocs_info <- arrow::read_feather(str_c(dir_extracted, "/", blockname_process, "_epocs_info.feather"))
    }

    print(" - extracting events and determining events of interest for peth")

    # events for perievent time histogram

    # get events_of_interest
    events_of_interest <- key_file$events_of_interest %>% str_split(";")
    events_of_interest <- events_of_interest[[1]] %>% str_trim()

    print(str_c(" - using events defined in log_fp: ", str_c(events_of_interest, collapse = ", ")))

    if (flag_events_tdt_epocs) {
      print("   ~ sourcing events from TDT epocs")

      events <-
        fp_epocs_to_events(
          epocs_data,
          key_file
        ) %>%
        group_by(event_id_char) %>%
        mutate(event_number = row_number()) %>%
        ungroup()


      if (nrow(events) == 0) {
        print("   * NO EVENTS DETECTED in tdt epoch data *")
        flag_events_tdt_epocs <- 0
      } else {
        ### return filtered events based on filters defined in fp_events_of_interest
        if (sum(!is.na(events_of_interest)) > 0) {
          # return empty filtered events if events of interest are listed in log_fp


          ## use clock to align arduino data to fp data ----------------------------------------------------------------
          if ("clock" %in% events$event_id_char %>% unique()) {
            print("    ~ sourcing events from arduino file synched to TDT clock")

            if (flag_multi_subject) {
              blockname_alignment <- blockname_process_multi_subject
            } else {
              blockname_alignment <- blockname_process
            }

            fn_arduino_events <- str_c(dir_extracted_arduino, blockname_alignment, "_event.csv")

            if (!file.exists(fn_arduino_events)) {
              print("!error: no extracted arduino file exists")
              print(fn_arduino_events)
              return(0)
            }

            df_beh_events_and_tick <- read.csv(fn_arduino_events) %>%
              mutate(event_ts = event_ts / 1000) %>%
              rename(ts = event_ts)

            ticks <- df_beh_events_and_tick %>%
              filter(event_id_char == "tick_clock_out") %>%
              mutate(tick_number = row_number())

            df_rec_system_tick <- events %>%
              filter(event_id_char == "clock") %>%
              tail(nrow(ticks)) %>%
              select(ts = event_ts) %>%
              mutate(event_id_char = "tick_clock_in")

            events <- convert_beh_arduino_to_fp_time(df_rec_system_tick, df_beh_events_and_tick)

            events <- events %>%
              filter(!is.na(ts_rec_system)) %>%
              mutate(blockname = blockname_process) %>%
              rename(event_ts = ts_rec_system, event_ts_arduino = ts) %>%
              group_by(event_id_char) %>%
              mutate(event_number = row_number()) %>%
              select(blockname, event_id_char, event_ts, event_ts_arduino, event_number)
          }



          events_filtered <- events %>%
            fp_filter_events(events_of_interest)

          flag_events_filtered <- 1
        } else {
          print("   * NO SPECIFICED EVENTS IN log_fp$events_of_interest *")
        }
      }
    } else {
      print("   * no TDT epocs file exists")
    }

    # If there is manual time stamps, read and join with filtered events
    events_manual_fn <- key_file$events_manual_fn # get events_manual_fn

    if (!is.na(events_manual_fn)) {
      print(str_c("   ~ sourcing evens from manual file: ", events_manual_fn))

      events_manual <- read_csv(events_manual_fn) %>%
        filter(blockname == blockname_process) %>%
        mutate(event_id_ttl = "manual")

      if (nrow(events_manual) == 0) {
        print("   * NO EVENTS DETECTED in specified manual file *")
      } else {
        events_manual_filtered <- events_manual %>%
          fp_filter_events(events_of_interest)

        if (flag_events_filtered) {
          events <- events %>% bind_rows(events_manual)
          events_filtered <- events_filtered %>% bind_rows(events_manual_filtered)
        } else {
          events <- events_manual
          events_filtered <- events_manual_filtered
          flag_events_filtered <- 1
        }

        if (nrow(events_manual_filtered) == 0) {
          print("   * manual events file found but no events made it through filter")
        }
      }
    } else {
      print("   * no specified manual events file in log_fp")
    }

    if (flag_events_filtered) {
      if (nrow(events_filtered) == 0) {
        flag_events_filtered <- 0
        print(str_c(" * no perievent time histogram: no filtered time stamps"))
      }
    }

    # perform perievent time histogram is there are events detected
    if (flag_events_filtered) {
      print(str_c(" - computing peth for ", dim(events_filtered)[1], " events"))

      events_filtered <- events_filtered %>%
        group_by(event_id_char) %>%
        mutate(event_number = row_number()) %>%
        ungroup()

      streams_data_smoothed_fitted_peth_downsampled <-
        fp_peri_event_time_histogram(
          streams_data_smoothed_fitted,
          events_filtered,
          key_file$peth_pre,
          key_file$peth_post,
          fs,
          key_file$down_sampled_freq
        )

      # }
    }

    # save processed data
    print(str_c(" - saving files to directory ", dir_processed))

    if (flag_multi_subject) {
      blockname_process <- blockname_process_multi_subject

      streams_data_smoothed_fitted_downsampled <- streams_data_smoothed_fitted_downsampled %>% mutate(blockname_multi_subject = blockname_process)

      events <- events %>% mutate(blockname_multi_subject = blockname_process)

      events_filtered <- events_filtered %>% mutate(blockname_multi_subject = blockname_process)

      streams_data_smoothed_fitted_peth_downsampled <- streams_data_smoothed_fitted_peth_downsampled %>% mutate(blockname_multi_subject = blockname_process)
    }

    print(str_c("   ~ ", blockname_process, "_streams_session.feather"))

    streams_data_smoothed_fitted_downsampled %>%
      mutate(blockname = blockname_process) %>%
      arrow::write_feather(str_c(dir_processed, "/", blockname_process, "_streams_session.feather"))

    if (flag_events_filtered) { # save events and peth if they exist
      print(str_c("   ~ ", blockname_process, "_events.feather"))

      events %>%
        mutate(blockname = blockname_process) %>%
        arrow::write_feather(str_c(dir_processed, "/", blockname_process, "_events.feather"))

      if (dim(events_filtered)[1] > 0) {
        print(str_c("   ~ ", blockname_process, "_events_peth.feather"))

        events_filtered %>%
          mutate(blockname = blockname_process) %>%
          arrow::write_feather(str_c(dir_processed, "/", blockname_process, "_events_peth.feather"))

        print(str_c("   ~ ", blockname_process, "_streams_peth.feather"))

        streams_data_smoothed_fitted_peth_downsampled %>%
          mutate(blockname = blockname_process) %>%
          arrow::write_feather(str_c(dir_processed, "/", blockname_process, "_streams_peth.feather"))
      }
    }
  }

  print("")
  print("fp_preprocess() complete")
}

# misc functions -------------------------------------------------------------------------------------------------------

return_channel_config <- function(streams) {
  # return the configuration of the recording session
  # - one fiber / one signal wavelength
  # - two fiber / one signal wavelength
  # - two fiber / two signal wavelengths
  #
  # this config will be used to determine which summary plot is made

  channel_info <- streams %>%
    select(fiber_id, signal_wavelength) %>%
    unique()

  num_fiber_id <- channel_info %>%
    select(fiber_id) %>%
    unique() %>%
    nrow()
  num_signal_wavelength <- channel_info %>%
    select(signal_wavelength) %>%
    unique() %>%
    nrow()

  if (num_fiber_id == 1 & num_signal_wavelength == 1) {
    channel_config <- "one_fiber_one_signal"
  }

  if (num_fiber_id == 2 & num_signal_wavelength == 1) {
    channel_config <- "two_fiber_one_signal"
  }

  if (num_fiber_id == 1 & num_signal_wavelength == 2) {
    channel_config <- "one_fiber_two_signal"
  }

  if (num_fiber_id == 2 & num_signal_wavelength == 2) {
    channel_config <- "two_fiber_two_signal"
  }

  return(channel_config)
}


filter_data_fp <- function(df, subjects_to_remove, subject_x_fibers_to_remove) {
  # filter out subjects and specific fibers

  # df (dataframe)
  #   - must contain subject_id for subject removal
  #   - must contain subject_id and channel_id for subject removal
  # subjects_to_remove (vector of strings), each subject to be removed (e.g. c('sbj01', 'sbj03'))
  # subject_x_fibers_to_remove (list of vectors of strings)
  #   - each list item must include a subject id and channel id to be removed (e.g. list(c('sbj01', 'channel_id01')))

  if (sum(str_detect(df %>% colnames(), "subject")) == 0) {
    print("no subject variable in df, cannot perform filter_data_fp()")
    return()
  }

  # remove subjects
  for (subject_id in subjects_to_remove) {
    df <- df %>%
      filter(subject != subject_id)
  }

  # remove specific fibers from a subject
  if (sum(str_detect(df %>% colnames(), "channel_id")) > 0) {
    if (sum(!is.na(subject_x_fibers_to_remove)) > 0) {
      for (n_pair in 1:length(subject_x_fibers_to_remove)) {
        df <- df %>%
          filter(!(subject == subject_x_fibers_to_remove[[n_pair]][1] &
            channel_id == subject_x_fibers_to_remove[[n_pair]][2])) # removed for poor signal in NAc lat
      }
    }
  }

  return(df)
}

fiber_id_to_channel_id <- function(df, log_fp) {
  # convert fiber_id from tdt data to channel_id as specififed in log_fp
  #
  # fiber_id's must be columns in log_fp starting with 'fiber_id' (e.g. fiber_id01, fiber_id02)

  df <- df %>%
    mutate(fiber_id = fiber_id %>% str_c("fiber_id0", .)) %>%
    left_join(log_fp %>%
      select(blockname, starts_with("fiber_id")) %>%
      gather(fiber_id, channel_id, starts_with("fiber_id")))

  return(df)
}

convert_blockname_multi_subject <- function(df) {
  # convert blockname_multi_subject to blockname for files with multi-subjects
  df %>%
    mutate(blockname = ifelse(is.na(blockname_multi_subject), blockname, blockname_multi_subject))
}


return_stim_params <- function(import_fns, dir_extracted, file_format_output) {
  # takes extracted serial data and returns stim parameters for individual stimulation events
  #
  # import_fns (vector of strings): vector of blocknames used for combined_import
  # dir_extracted (string): path to the directory for serial extraction
  #
  # output: stim_params
  #   1 row / stimulation with columns for each parameter (defined in filter below)

  data_list <- combined_import(
    import_directory = dir_extracted,
    filter_strings = import_fns,
    prefixes = NA,
    suffixes = c(str_c("_param_dynamic.", file_format_output))
  )

  stim_params <- data_list[[1]] %>%
    select(blockname, param_dynamic, param_ts, param_value) %>%
    filter(param_dynamic %in% c("trial_stim_frequency", "trial_stim_pulse_duration", "trial_stim_train_duration")) %>%
    arrange(blockname, param_ts) %>%
    spread(param_dynamic, param_value) %>%
    group_by(blockname) %>%
    mutate(stim_num = row_number())

  return(stim_params)
}


return_local_zscore <- function(streams_peth, time_window) {
  # computes zscore relative to mean and standard deviation during time_window prior to events
  #
  # inputs
  #  - streams_peth (dataframe): must contain blockname, event_id_char, time_rel, and delta_signal_poly
  #  - time_window (vector, double): time in sec. for window used for mu / sd (e.g. c(-3, 0) for 3s baseline )
  #
  # outputs:
  #  - streams_peth with new columns
  #     ~ delta_signal_poly_zscore_local
  #     ~ delta_control_poly_zscore_local

  streams_peth <- streams_peth %>%
    group_by(blockname, event_id_char) %>%
    filter(time_rel > time_window[1], time_rel < time_window[2]) %>%
    summarise(
      delta_signal_poly_local_bl_mu = mean(delta_signal_poly),
      delta_signal_poly_local_bl_sd = sd(delta_signal_poly),
      delta_control_poly_local_bl_mu = mean(delta_control_poly),
      delta_control_poly_local_bl_sd = sd(delta_control_poly),
      .groups = "drop"
    ) %>%
    left_join(streams_peth, ., by = c("blockname", "event_id_char")) %>%
    mutate(
      delta_signal_poly_zscore_local = (delta_signal_poly - delta_signal_poly_local_bl_mu) / delta_signal_poly_local_bl_sd,
      delta_control_poly_zscore_local = (delta_control_poly - delta_control_poly_local_bl_mu) / delta_control_poly_local_bl_sd
    ) %>%
    select(-delta_signal_poly_local_bl_mu, -delta_signal_poly_local_bl_sd)

  return(streams_peth)
}

inspect_raw_fp_folder <- function(dir_raw) {
  dirs <- list.dirs(dir_raw)

  dirs <- dirs[2:length(dirs)]

  print("checking for subfolders and folders with incorrect number of files in raw fp folder...")
  print("")
  # determine number of folders with subfolders within raw fp directory
  n_subfolders


  dirs_subfolders <- dirs %>%
    as_tibble() %>%
    mutate(value = str_remove(value, dir_raw)) %>%
    rowwise() %>%
    mutate(n_subfolders = str_count(value, "/") - 1) %>%
    filter(n_subfolders > 0)


  if (nrow(dirs_subfolders) == 0) {
    print("no erroneous subfolders in folders within raw directory")
  } else {
    print("erroneous subfolders...")
    for (folder in 1:nrow(dirs_subfolders)) {
      print(dirs_subfolders$value[folder])
    }
  }
  print("")

  # determine number of files within each folder within raw fp directory
  n_file <- c()

  for (d in dirs) {
    n_file <- length(list.files(d)) %>% c(n_file, .)
  }

  dirs_file_counts <- dirs %>%
    as_tibble() %>%
    mutate(n = n_file) %>%
    filter(n != 8)

  if (nrow(dirs_file_counts) == 0) {
    print("no folders with incorrect number of files")
  } else {
    print("folders with incorrect number of files...")
    for (folder in 1:nrow(dirs_file_counts)) {
      print(dirs_file_counts$value[folder])
    }
  }
}


# scaling external events ----------------------------------------------------------------------------------------------
return_tdt_external_delta <- function(events_tdt, events_external, key_var) {
  tdt <- events_tdt %>%
    ungroup() %>%
    filter(event_id_char %in% key_var) %>%
    select(blockname, event_id_char, event_ts) %>%
    group_by(blockname) %>%
    mutate(trim_time_tdt = min(event_ts)) %>%
    mutate(event_ts = event_ts - trim_time_tdt) %>%
    arrange(blockname, event_ts, event_id_char) %>%
    mutate(event_number = row_number()) %>%
    rename(event_ts_tdt = event_ts)

  if (nrow(tdt) == 0) {
    print(str_c("no events in tdt that match key_var: ", key_var))
  }

  external <- events_external %>%
    filter(event_id_char %in% c(key_var)) %>%
    filter(event_ts > 0) %>%
    select(blockname, event_id_char, event_ts) %>%
    mutate(event_ts = event_ts) %>%
    group_by(blockname) %>%
    mutate(trim_time_external = min(event_ts)) %>%
    mutate(event_ts = event_ts - trim_time_external) %>%
    arrange(blockname, event_ts, event_id_char) %>%
    rename(event_ts_external = event_ts)

  if (nrow(external) == 0) {
    print(str_c("no events in external that match key_var: ", key_var))
  }

  if (nrow(tdt) == 0 | nrow(external) == 0) {
    return()
  }

  print(left_join(
    tdt %>%
      group_by(blockname) %>%
      summarise(n_row_tdt = n()),
    external %>%
      group_by(blockname) %>%
      summarise(n_row_external = n())
  ))

  scale_factors <- tdt %>%
    left_join(external) %>%
    group_by(blockname) %>%
    mutate(step_factor = event_ts_tdt - event_ts_external) %>%
    mutate(event_ts = event_ts_external + trim_time_external) %>%
    select(blockname, event_id_char, event_ts, trim_time_tdt, trim_time_external, step_factor)

  return(scale_factors)
}

scale_external_to_tdt_local_delta <- function(events_tdt, events_external, key_var) {
  # convert external time to time in tdt file based on a set of shared events
  # inputs. Uses a simple delta to shift external events to match timing of tdt events.
  # by filling down, the delta of key events can be used to shift following events not
  # found in the tdt file.

  # - events_tdt (dataframe) events from external (time unit needs to be consistent with tdt (sec))
  # - events_external (dataframe) events from tdt
  # - key_var (vector of strings) strings of variable names to use to create local deltas
  #
  # variables within events_tdt and events_external
  #  - blockname (unique file identifier)
  #  - event_id_char (string id for event type, must match for 2 dataframes)
  #  - event_ts (time of event, units must match for 2 dataframes)

  scale_factors <- return_tdt_external_delta(events_tdt, events_external, key_var)

  data_rotary_adj <- events_external %>%
    left_join(scale_factors %>% select(blockname, trim_time_tdt, trim_time_external) %>% unique()) %>%
    left_join(scale_factors %>% select(blockname, event_id_char, event_ts, step_factor) %>% unique()) %>%
    arrange(blockname, event_ts) %>%
    group_by(blockname) %>%
    fill(trim_time_tdt, trim_time_external, step_factor) %>% # fill down steps and trims to adjust data not recorded by tdt
    mutate(event_ts_tdt = event_ts - trim_time_external + trim_time_tdt) %>% # shift time to be relative to tdt onset
    mutate(event_ts_tdt = event_ts_tdt + step_factor) # shift data using local delta

  return(data_rotary_adj)
}


return_tdt_external_scaled <- function(events_tdt, events_external, key_vars) {
  tdt <- events_tdt %>%
    ungroup() %>%
    filter(event_id_char %in% key_vars) %>%
    select(blockname, event_id_char, event_ts) %>%
    group_by(blockname) %>%
    mutate(trim_time_tdt = min(event_ts)) %>%
    mutate(event_ts = event_ts - trim_time_tdt) %>%
    filter(event_ts == max(event_ts)) %>%
    unique() %>%
    rename(event_ts_tdt = event_ts) %>%
    arrange(blockname, event_ts_tdt) %>%
    group_by(blockname)

  if (nrow(tdt) == 0) {
    print(str_c("no events in tdt that match key_vars: ", key_vars))
  }

  external <- events_external %>%
    filter(event_id_char %in% c(key_vars)) %>%
    filter(event_ts > 0) %>%
    select(blockname, event_id_char, event_ts) %>%
    mutate(event_ts = event_ts) %>%
    group_by(blockname) %>%
    mutate(trim_time_external = min(event_ts)) %>%
    mutate(event_ts = event_ts - trim_time_external) %>%
    filter(event_ts == max(event_ts)) %>%
    unique() %>%
    rename(event_ts_external = event_ts)

  if (nrow(external) == 0) {
    print(str_c("no events in external that match key_vars: ", key_vars))
  }

  if (nrow(tdt) == 0 | nrow(external) == 0) {
    return()
  }

  scale_factors <- tdt %>%
    left_join(external) %>%
    mutate(scale_factor = (event_ts_external - event_ts_tdt) / event_ts_tdt) %>%
    select(blockname, trim_time_external, trim_time_tdt, scale_factor)

  return(scale_factors)
}

scale_external_to_tdt_session_scale <- function(events_tdt, events_external, key_vars) {
  # convert external time to time in tdt file by scaling external time based on the
  # first and last events defined by key_vars.
  #
  # (total time external - total time tdt) / total time tdt)

  # - events_tdt (dataframe) events from external (time unit needs to be consistent with tdt (sec))
  # - events_external (dataframe) events from tdt
  # - key_var (vector of strings) strings of variable names to use to create global scale factor
  #
  # variables within events_tdt and events_external
  #  - blockname (unique file identifier)
  #  - event_id_char (string id for event type, must match for 2 dataframes)
  #  - event_ts (time of event, units must match for 2 dataframes)

  scale_factors <- return_tdt_external_scaled(events_tdt, events_external, key_vars)

  data_rotary_adj <- events_external %>%
    left_join(scale_factors) %>%
    arrange(blockname, event_ts) %>%
    mutate(event_ts_tdt = event_ts - trim_time_external) %>% # shift time to be relative to first event in external
    mutate(event_ts_tdt = event_ts_tdt - event_ts_tdt * scale_factor) %>% # scale time based on session scale factor
    mutate(event_ts_tdt = event_ts_tdt + trim_time_tdt) # shift time to tdt start

  return(data_rotary_adj)
}

# resampling time series data ------------------------------------------------------------------------------------------

get_response_set_rotary_intermittent <- function(df, event_id_char_setend) {
  # returns events for break_engaged following active rotation, and active_rotation_criteria
  # and generates a new variable that defines response sets of breaks and the subsequent criteria
  #  (break, ..., active_rotation_criteria)
  #
  # used for resampling time series data inbetween each individual response

  event_set <- df %>%
    arrange(blockname, event_ts, event_id_char) %>%
    unique() %>%
    filter(!is.na(event_ts_tdt)) %>%
    group_by(blockname, event_id_char) %>%
    mutate(response_set = row_number()) %>%
    mutate(response_set = ifelse(event_id_char == event_id_char_setend, response_set, NA)) %>%
    group_by(blockname) %>%
    fill(response_set, .direction = "up") %>%
    filter(!is.na(response_set)) # filter out responses that are not followed by criteria

  return(event_set)
}

resample_time_series_fp <- function(streams, event_set, tm_pre, tm_post, resample_rate, ...) {
  # return resampled time series data relative to sets of events using resample_rate
  # the time between each event, time prior to first event, and time following last event
  # will be each be set to 1, with the total time equal to 1 + the number of events

  # streams: fp stream data
  # event_set: set of events to define resampling stakes
  #   - must contain event_ts_tdt for timeing of events relative to tdt
  #   - must contain response_set for defining each set of events
  # tm_pre: time prior to first event
  # tm_post: time following last event
  # resample_rate: samping rate within each resampled period within the data

  blocknames <- event_set %>% # get blocknames from event_set
    select(blockname) %>%
    unique() %>%
    pull()

  first_loop <- 1

  for (blockname_loop in blocknames) { # for each blockname
    print(blockname_loop)

    rotation_event_sets_loop <- event_set %>% # filter event set for blockname
      filter(blockname == blockname_loop) %>%
      arrange(event_ts_tdt)

    rescale_streams_loop <- streams %>% # filter streams for blockname
      filter(blockname == blockname_loop)

    for (n_set in seq(1, rotation_event_sets_loop %>% pull(response_set) %>% max(na.rm = TRUE))) { # for each event set

      rotation_event_set <- rotation_event_sets_loop %>% # generate flanks of resamping period
        filter(response_set == n_set) %>%
        select(blockname, event_ts_tdt) %>%
        mutate(event_ts_tdt_end = lead(event_ts_tdt))


      # for first period retrieve data prior to first event and convert it to proportion of time period
      rescale_streams_resampled <- rescale_streams_loop %>%
        filter(time > rotation_event_set$event_ts_tdt[1] - tm_pre, time < rotation_event_set$event_ts_tdt[1]) %>%
        mutate(time_resampled = time - min(time)) %>%
        mutate(time_resampled = time_resampled / max(time_resampled)) %>%
        mutate(resample_period = 0)

      # for periods between flanks, retrieve data between flanks and convert it to proportion of time period
      for (n_period in 1:(nrow(rotation_event_set) - 1)) {
        rescale_streams_resampled <- rescale_streams_loop %>%
          filter(time >= rotation_event_set$event_ts_tdt[n_period], time < rotation_event_set$event_ts_tdt_end[n_period]) %>%
          mutate(time_resampled = time - min(time)) %>%
          mutate(time_resampled = time_resampled / max(time_resampled)) %>%
          mutate(resample_period = n_period) %>%
          bind_rows(rescale_streams_resampled, .)
      }

      # for last period retrieve data following last event and convert it to proportion of time period
      rescale_streams_resampled <- rescale_streams_loop %>%
        filter(time >= rotation_event_set$event_ts_tdt_end[n_period], time < rotation_event_set$event_ts_tdt_end[n_period] + tm_post) %>%
        mutate(time_resampled = time - min(time)) %>%
        mutate(time_resampled = time_resampled / max(time_resampled)) %>%
        mutate(resample_period = n_period + 1) %>%
        bind_rows(rescale_streams_resampled, .)

      # use resample_period to shift data into a continuous stream
      rescale_streams_resampled <- rescale_streams_resampled %>%
        mutate(time_resampled = time_resampled + resample_period)

      # generate time bins
      tm_bins <- seq(
        -1,
        rescale_streams_resampled %>% pull(time_resampled) %>% max(na.rm = TRUE),
        resample_rate
      )

      # resample data to generate average activity during each of the newly defined samples proportional to events
      rescale_streams_resampled <- rescale_streams_resampled %>%
        mutate(time_rel = cut(time_resampled, tm_bins, tm_bins[2:length(tm_bins)])) %>%
        select(blockname, fiber_id, signal_wavelength, control_wavelength, time_rel, time:delta_control_poly_dff) %>% #*** variables defined here ***
        group_by(blockname, fiber_id, signal_wavelength, control_wavelength, time_rel) %>%
        summarise_all(mean, .groups = "drop") %>% # summarise all variables using mean
        mutate(time_rel = time_rel %>% as.character() %>% as.double()) %>% # convert time to numeric
        filter(time_rel < max(time_rel)) %>%
        mutate(n_response_set = n_set)

      # combine data from each set / blockname into a single dataframe
      if (first_loop) {
        rescale_streams_resampled_out <- rescale_streams_resampled
        first_loop <- 0
      } else {
        rescale_streams_resampled_out <- rescale_streams_resampled %>% bind_rows(rescale_streams_resampled_out, .)
      }
    }
  }
  return(rescale_streams_resampled_out) # return combined dataframe
}

# neurophotometrics preprocessing --------------------------------------------------------------------------------------
npm_extract <- function(log_npm, dir_npm_raw, dir_npm_extracted, manual_session_names, overwrite, filter_red_opto_sd) {
  if (missing(filter_red_opto_sd)) {
    filter_red_opto_sd <- NA
  }


  # sessions in dir_npm_raw
  raw_sessions <- list.files(dir_npm_raw)

  # generate sessions_to_process
  sessions_to_process <- tibble("session_name" = raw_sessions)


  # file name suffixes in dir_npm_extracted
  suffixes <- c(
    "_clock.feather",
    "_streams_data.feather",
    "_session_info.feather",
    "_plt_trims.png"
  )

  if (overwrite == 0) {
    fns_diff <- dir_diff(dir_npm_raw, dir_npm_extracted, suffixes, print_message = 0)

    sessions_to_process <- sessions_to_process %>%
      filter(session_name %in% fns_diff)

    # filter block names produced by multi_subject

    if (nrow(sessions_to_process) == 0) {
      print("No new data to process")
      print("")
      print("to overwrite data, set overwrite = 1")

      return()
    }
  } else {
    print("overwriting data in dir_npm_extracted...")
    print("")
  }

  # filter if using manual_session_names
  if (sum(!is.na(manual_session_names)) > 0) {
    sessions_to_process <- sessions_to_process %>%
      filter(session_name %in% manual_session_names)

    if (nrow(sessions_to_process) == 0) {
      print("no data in dir_raw matching manual_session_names")
      print("")

      return()
    }
  }

  print("processing following session_name:")
  for (session_name in sessions_to_process$session_name) {
    print(session_name)
  }
  print("")

  # loop over each session_name_________________________________________________________________________________
  for (session_name_process in sessions_to_process$session_name) { # loop over each session
    print(str_c("extracting session_name: ", session_name_process))

    # import data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    data_npm_fp <- read.csv(str_c(dir_npm_raw, session_name_process, "/", "data_fp.csv")) %>%
      clean_names() %>%
      mutate(session_name = session_name_process)

    data_npm_fp_ts <- read.csv(str_c(dir_npm_raw, session_name_process, "/", "data_fp_ts.csv")) %>%
      clean_names() %>%
      mutate(session_name = session_name_process)

    data_npm_serial <- read.csv(str_c(dir_npm_raw, session_name_process, "/", "data_serial.csv")) %>%
      clean_names() %>%
      mutate(session_name = session_name_process)


    # filter frames with red opto interference (optional)
    if (!is.na(filter_red_opto_sd)) {
      data_npm_fp <- remove_opto_stim_frames_npm(data_npm_fp, filter_red_opto_sd)
    }

    # join ts_comp  from data_npm_fp_ts to data_npm_fp
    data_npm_fp <- data_npm_fp %>%
      left_join(
        data_npm_fp_ts %>%
          rename(
            frame_counter = item1,
            ts_comp = item2
          ),
        by = c("frame_counter", "session_name")
      )


    # channel info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    var_channels <- c("channel_wavelength", "channel_power", "led_state", "channel_color")
    var_join <- c("session_name", "blockname", "branch_id", "region", "down_sampled_freq", "channel_number")

    log_npm_filtered_to_channels <- log_npm %>%
      select(session_name, blockname, branch_id, region, down_sampled_freq, all_of(var_channels))

    for (var_channel in var_channels) {
      loop_df <- log_npm_filtered_to_channels %>%
        select(session_name, blockname, branch_id, region, down_sampled_freq, !!as.name(var_channel)) %>%
        separate(!!(as.name(var_channel)), into = c("ch01", "ch02", "ch03"), fill = "right") %>%
        gather("channel_number", !!var_channel, ch01:ch03)

      if (var_channel == var_channels[1]) {
        log_npm_long <- loop_df
      } else {
        log_npm_long <- loop_df %>%
          left_join(log_npm_long, ., by = var_join)
      }
    }

    log_npm_long <- log_npm_long %>%
      mutate(roi = ifelse(channel_color == "r", (branch_id - 1) * 2, (branch_id * 2) - 1)) %>%
      mutate(led_state = led_state %>% as.integer()) %>%
      filter(session_name == session_name_process)


    data_npm_fp_long <- data_npm_fp %>%
      gather("roi", "raw_au", starts_with("region")) %>%
      mutate(roi = roi %>% str_remove("region")) %>%
      mutate(
        channel_color = roi %>% str_sub(-1),
        roi = roi %>% str_remove(channel_color) %>% as.integer()
      )

    # filter to data_npm_fp_long to signals defined in log_npm_long
    data_npm_fp_long_filtered <- data_npm_fp_long %>%
      inner_join(log_npm_long, by = c("session_name", "led_state", "roi", "channel_color"))


    # filter timeseries based on serial clock input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    flag_clock <- log_npm %>%
      filter(session_name == session_name_process) %>%
      pull(clock_pin)

    flag_clock <- sum(!is.na(flag_clock)) > 0

    if (flag_clock) { # if there is a clock noted in log_npm
      data_npm_serial_clock <- data_npm_serial %>%
        arrange(ts_comp) %>%
        unique() %>%
        group_by(clock_pin) %>%
        mutate(tick_number = row_number()) %>%
        ungroup() %>%
        left_join(log_npm %>% select(session_name, blockname, clock_pin, tm_baseline, tm_post) %>% unique(),
          by = c("clock_pin", "session_name")
        ) %>%
        filter(!is.na(blockname))
    }

    # pull out first and last clock
    if (flag_clock) {
      filter_times <- data_npm_serial_clock %>%
        group_by(session_name, blockname, clock_pin) %>%
        filter(ts_comp == min(ts_comp) | ts_comp == max(ts_comp)) %>%
        mutate(break_id = ifelse(ts_comp == min(ts_comp), "tm_session_start", "tm_session_end")) %>%
        mutate(ts_comp = ifelse(break_id == "tm_session_start",
          ts_comp - tm_baseline * 1000,
          ts_comp + tm_post * 1000
        )) %>%
        ungroup() %>%
        unique()
    } else {
      filter_times <- data_npm_fp_long_filtered %>%
        left_join(log_npm %>% select(session_name, blockname, tm_trim_start, tm_trim_end)) %>%
        group_by(session_name, blockname) %>%
        filter(ts_comp > min(ts_comp) + (tm_trim_start * 1000) | ts_comp < max(ts_comp) - tm_trim_end * 1000) %>%
        filter(ts_comp == min(ts_comp) | ts_comp == max(ts_comp)) %>%
        mutate(break_id = ifelse(ts_comp == min(ts_comp), "tm_session_start", "tm_session_end")) %>%
        ungroup() %>%
        select(ts_comp, session_name, blockname, break_id) %>%
        unique()
    }

    # filter data based on first and last clock
    data_npm_fp_long_filtered <- data_npm_fp_long_filtered %>%
      left_join(
        # wide breaks for session start and end
        filter_times %>%
          select(session_name, blockname, break_id, ts_comp) %>%
          unique() %>%
          spread(break_id, ts_comp),
        by = c("session_name", "blockname")
      ) %>%
      filter(ts_comp > tm_session_start, ts_comp < tm_session_end)

    # resample and interpolate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    data_npm_fp_long_filtered_resampled <- data_npm_fp_long_filtered %>%
      select(session_name, blockname, branch_id, channel_wavelength, down_sampled_freq, ts_comp, raw_au) %>%
      group_by(session_name, blockname, branch_id, channel_wavelength) %>%
      group_modify(~ timeseries_resample(., "ts_comp", "raw_au")) %>%
      group_modify(~ timeseries_interpolate_rows(., "ts_comp", "raw_au")) %>%
      ungroup()

    # find event times ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    data_npm_input_events <- data_npm_fp_long_filtered %>%
      select(session_name, blockname, ts_comp, input0, input1) %>%
      gather(input_id, state, input0:input1) %>%
      arrange(input_id, ts_comp) %>%
      group_by(input_id) %>%
      mutate(event = ifelse(state == 1 & lag(state == 0), "event_ts", NA)) %>%
      mutate(event = ifelse(state == 0 & lag(state == 1), "event_ts_offset", event)) %>%
      filter(!is.na(event)) %>%
      group_by(session_name, blockname, input_id, event) %>%
      mutate(event_number = row_number()) %>%
      select(-state) %>%
      spread(event, ts_comp)


    # generate npm / computer timestamp dataframe ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    time_stamp_key <- data_npm_fp %>%
      select(session_name, frame_counter, timestamp, ts_comp) %>%
      unique() %>%
      tibble() %>%
      rename(
        ts_npm = timestamp,
        ts_comp = ts_comp
      )

    # generate trim plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ts_comp_subtract <- data_npm_fp_long_filtered_resampled$ts_comp %>% min()

    plt_trims <- data_npm_fp_long_filtered_resampled %>%
      slice_sample(prop = 0.01) %>%
      mutate(ts_comp = (ts_comp - ts_comp_subtract) / 1000 / 60) %>%
      ggplot(aes(ts_comp, raw_au, color = channel_wavelength, group = interaction(blockname, channel_wavelength, branch_id))) +
      geom_line() +
      geom_vline(data = filter_times, aes(xintercept = (ts_comp - ts_comp_subtract) / 1000 / 60), linetype = 2) +
      facet_grid(blockname ~ ., scales = "free") +
      theme_minimal() +
      xlab("time (m)") +
      ylab("raw (au)") +
      force_panelsizes(cols = unit(3.5, "cm"), rows = unit(3.5, "cm"))

    # generate on plot
    if (flag_clock) {
      last_tm_session_start <- filter_times %>%
        ungroup() %>%
        filter(break_id == "tm_session_start") %>%
        filter(ts_comp == max(ts_comp)) %>%
        pull(ts_comp)

      plt_starts <- data_npm_fp_long_filtered_resampled %>%
        slice_sample(prop = 0.1) %>%
        filter(ts_comp < last_tm_session_start + 10 * 1000) %>%
        mutate(ts_comp = (ts_comp - ts_comp_subtract) / 1000 / 60) %>%
        group_by(blockname, branch_id, ts_comp) %>%
        summarise(raw_au = mean(raw_au), .groups = "drop") %>%
        ggplot(aes(ts_comp, raw_au, color = blockname, group = interaction(blockname, branch_id))) +
        geom_line() +
        geom_vline(
          data = filter_times %>% filter(break_id == "tm_session_start"),
          aes(
            xintercept = (ts_comp - ts_comp_subtract) / 1000 / 60,
            linetype = blockname
          )
        ) +
        theme_minimal() +
        xlab("time (m)") +
        ylab("mean raw (au)") +
        force_panelsizes(cols = unit(3.5, "cm"), rows = unit(3.5, "cm"))

      ggsave(str_c(dir_npm_extracted, session_name_process, "_plt_starts.png"), plt_starts, width = unit(5, "cm"), height = unit(5, "cm"))
    }

    # export ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    log_npm_long %>%
      left_join(log_npm %>% select(blockname, clock_pin)) %>%
      write_feather(str_c(dir_npm_extracted, session_name_process, "_session_info.feather"))

    data_npm_fp_long_filtered_resampled %>%
      write_feather(str_c(dir_npm_extracted, session_name_process, "_streams_data.feather"))

    if (flag_clock) {
      data_npm_serial_clock %>%
        select(session_name, blockname, clock_pin, tick_number, ts_clock, ts_comp) %>%
        write_feather(str_c(dir_npm_extracted, session_name_process, "_clock.feather"))
    }

    if (nrow(data_npm_input_events) > 0) {
      data_npm_input_events %>%
        write_feather(str_c(dir_npm_extracted, session_name_process, "_input_data.feather"))
    }

    time_stamp_key %>%
      write_feather(str_c(dir_npm_extracted, session_name_process, "_timestamps.feather"))


    ggsave(str_c(dir_npm_extracted, session_name_process, "_plt_trims.png"), plt_trims, width = unit(8, "cm"), height = unit(24, "cm"))


    print(log_npm_long %>%
      left_join(log_npm %>% select(blockname, clock_pin)) %>%
      select(session_name, blockname, clock_pin, branch_id) %>%
      unique() %>%
      group_by(session_name, blockname, clock_pin) %>%
      summarise(branch_id = paste0(branch_id, collapse = ";"), .groups = "drop") %>%
      left_join(filter_times %>%
        select(session_name, blockname, break_id, ts_comp) %>%
        spread(break_id, ts_comp)))
  } # loop over each session
}

fp_preprocess_npm <-
  function(dir_npm_extracted,
           dir_npm_processed,
           log_npm,
           manual_session_names = NA,
           overwrite = 0,
           dir_arduino_extracted = NA) {
    if (missing(dir_arduino_extracted)) {
      dir_arduino_extracted <- NA
    }


    extracted_files <- list.files(dir_npm_extracted)

    # file name suffixes in dir_npm_processed and dir_npm_extracted
    suffixes <- c(
      "_input_data.feather",
      "_streams_peth.feather",
      "_streams_session.feather",
      "_events_peth.feather",
      "_events.feather",
      "_clock.feather",
      "_streams_data.feather",
      "_session_info.feather",
      "_plt_trims.png",
      "_timestamps.feather",
      "_plt_starts.png"
    )

    # generate key
    key_files <- tibble("session_name" = extracted_files)

    if (nrow(key_files) == 0) {
      print("WARNING: No files in dir_npm_extracted (correct pathway?)")
      return()
    }

    # remove suffixes
    for (suffix in suffixes) {
      key_files <-
        key_files %>% mutate(session_name = session_name %>% str_remove(suffix))
    }

    key_files <- key_files %>% unique()

    # join info from log_npm
    key_files <- key_files %>%
      left_join(log_npm, by = "session_name")

    # filter if using manual_session_names
    if (sum(!is.na(manual_session_names)) > 0) {
      key_files <- key_files %>%
        filter(session_name %in% manual_session_names)

      if (nrow(key_files) == 0) {
        print("No data in dir_npm_extracted matching manual_session_names")
        print("")

        return()
      }
    }


    # check to determine if data in dir_npm_extracted are already located in dir_npm_processed
    if (overwrite == 0) {
      fns_diff <-
        dir_diff_session(
          dir_npm_extracted,
          dir_npm_processed,
          suffixes,
          print_message = 0,
          log_npm
        )

      key_files <- key_files %>%
        filter(session_name %in% fns_diff)

      if (nrow(key_files) == 0) {
        print("No new data to process")
        print("")
        print("to overwrite data, set overwrite = 1")

        return()
      }
    }


    sessions_to_process <- key_files %>%
      pull(session_name) %>%
      unique()


    # Bulk processing...
    print(str_c(
      "bulk processing ",
      length(sessions_to_process),
      " sessions..."
    ))

    for (session_process in sessions_to_process) {
      key_file <- key_files %>%
        filter(session_name == session_process)

      blocknames <- log_npm %>%
        filter(session_name == session_process) %>%
        pull(blockname) %>%
        unique()

      print(
        "_______________________________________________________________________________"
      )
      print("")
      print(str_c("processing session name: ", session_process))



      # determine signals within blockname
      signal_ids <- key_file %>%
        select(session_name, blockname, starts_with("signal_id")) %>%
        gather(
          "signal_id",
          "channel_wavelength",
          starts_with("signal_id")
        ) %>%
        filter(!is.na(channel_wavelength)) %>%
        arrange(signal_id) %>%
        left_join(key_file %>% select(session_name, blockname, control_id), by = c("session_name", "blockname"))


      # import extracted  data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      session_info <-
        arrow::read_feather(str_c(
          dir_npm_extracted,
          "/",
          session_process,
          "_session_info.feather"
        ))

      frame_time_stamps <-
        arrow::read_feather(str_c(
          dir_npm_extracted,
          "/",
          session_process,
          "_timestamps.feather"
        )) %>%
        mutate(ts_npm = ts_npm - min(ts_npm)) %>%
        mutate(ts_comp = ts_comp / 1000)

      # display session contents
      session_info_display <- session_info %>%
        select(blockname, branch_id) %>%
        unique() %>%
        group_by(blockname) %>%
        summarise(branch_id = paste0(branch_id, collapse = ";"))


      print("")
      print("- session contents:")

      for (nblock in seq(1, nrow(session_info_display))) {
        print(str_c(" - ", session_info_display$blockname[nblock], " - branch(s): ", session_info_display$branch_id[nblock]))
      }

      print("")

      print("- event sources:")

      # read in clock data if it exists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (file.exists(str_c(dir_npm_extracted, "/", session_process, "_clock.feather"))) {
        flag_clock <- 1
      } else {
        flag_clock <- 0
        print("   * no clock data")
      }

      if (flag_clock) {
        session_clock <-
          arrow::read_feather(str_c(
            dir_npm_extracted,
            "/",
            session_process,
            "_clock.feather"
          )) %>%
          mutate(ts_comp = ts_comp / 1000)
      }

      streams_data <-
        arrow::read_feather(str_c(
          dir_npm_extracted,
          "/",
          session_process,
          "_streams_data.feather"
        )) %>%
        mutate(ts_comp = ts_comp / 1000)


      streams_data <- streams_data %>%
        rename(time = ts_comp)

      # smooth data
      streams_data_smoothed <-
        fp_moving_average_grouped(streams_data, log_npm)


      # perform lls & poly corrections, zscores, and delta f/f
      streams_data_smoothed_fitted <-
        fp_normalization_npm(streams_data_smoothed, log_npm) %>%
        rename(signal_wavelength = channel_wavelength)




      # extract events-----------------------------------------------------------------------------------------------
      events <- tibble(
        session_name = character(),
        blockname = character(),
        input_id = character(),
        event_id_char = character(),
        event_ts = double(),
        event_ts_offset = double()
      )


      #  ~ npm inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      flag_events_npm_input <- 0

      # read in input data if it exists
      fn_npm_inputs <-
        str_c(
          dir_npm_extracted,
          "/",
          session_process,
          "_input_data.feather"
        )

      if (!file.exists(fn_npm_inputs)) {
        print("   * no npm input file")
        flag_events_npm_input <- 0
      } else {
        print("   ~ sourcing events from npm input")
        flag_events_npm_input <- 1
        npm_input_data <- arrow::read_feather(fn_npm_inputs)

        events_npm_input <-
          fp_npm_inputs_to_events(npm_input_data, log_npm) # convert npm_input_data to events_npm_input dataframe


        # combine events_npm_input with events dataframe
        events <- events %>%
          bind_rows(events_npm_input)
      }

      #  ~ manual events ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # If there is manual time stamps, read and join with filtered events
      events_manual_fns <-
        log_npm %>%
        filter(session_name == session_process) %>%
        filter(!is.na(events_manual_fn)) %>%
        pull(events_manual_fn) %>%
        unique()

      if (length(events_manual_fns) == 0) {
        print("   * no specified manual events file in log_npm")
      } else {
        print(str_c(
          "   ~ sourcing events from ",
          length(events_manual_fns),
          " manual events file(s)"
        ))

        # read in each specified manual events file for session and combine into one dataframe
        for (events_manual_fn in events_manual_fns) {
          # display warning and return if manual event file does not exist
          if (!file.exists(events_manual_fn)) {
            print("!error: events_manual_fn specified in log_npm does not exist")
            print(events_manual_fn)
            return(0)
          }

          # read in file and combine
          if (events_manual_fn == events_manual_fns[[1]]) {
            events_manual <- read.csv(events_manual_fn) %>%
              filter(session_name == session_process)
          } else {
            events_manual <- read.csv(events_manual_fn) %>%
              filter(session_name == session_process) %>%
              bind_rows(events_manual, .)
          }
        }

        # check to ensure that manual time stamps fall within recording time
        if (max(events_manual$event_ts) > max(frame_time_stamps$ts_npm)) {
          print("     * WARNING: manual event timestamp(s) exceed npm recording duration!")
        }

        # filter out time stamps that exceed recording duration
        events_manual <- events_manual %>%
          filter(event_ts < max(frame_time_stamps$ts_npm))

        # convert npm_time to ts_comp baed on frame_time_stamps
        events_manual <- bind_rows(
          events_manual %>%
            rename(ts_npm = event_ts),
          frame_time_stamps
        ) %>%
          arrange(ts_npm) %>%
          fill(frame_counter, .direction = "up") %>%
          fill(ts_comp, .direction = "up") %>%
          filter(!is.na(blockname)) %>%
          rename(event_ts = ts_comp) %>%
          select(-frame_counter, -ts_npm) %>%
          mutate(input_id = "manual")

        # combine events_manual with events dataframe
        events <- events %>%
          bind_rows(events_manual)
      }

      #  ~ clock events ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # If there is manual time stamps, read and join with filtered events
      key_clock <- log_npm %>%
        filter(session_name == session_process) %>%
        select(session_name, blockname, clock_pin) %>%
        unique()

      # if clocks are defined in log_npm
      if (sum(!is.na(key_clock$clock_pin)) == 0) {
        print("   * no clock inputs specified in log_npm")
      } else {
        print(str_c(
          "   ~ sourcing and aligning events from ",
          nrow(key_clock),
          " arduino file(s)"
        ))

        # read in each arduino file
        for (arduino_blockname in key_clock$blockname) {
          fn_arduino_events <-
            str_c(
              dir_arduino_extracted,
              arduino_blockname,
              "_event.csv"
            )

          # display warning and return if manual event file does not exist
          if (!file.exists(fn_arduino_events)) {
            print("!error: corresponding arduino behavior file does not exist")
            print(fn_arduino_events)
            return(0)
          }

          # read in data and combine
          if (arduino_blockname == key_clock$blockname[[1]]) {
            events_arduino <- read.csv(fn_arduino_events)
          } else {
            events_arduino <- read.csv(fn_arduino_events) %>%
              bind_rows(events_arduino, .)
          }
        }

        # save session_name
        events_arduino <- events_arduino %>%
          mutate(session_name = session_process) %>%
          select(session_name, everything())

        # convert time
        events_arduino <- events_arduino %>%
          mutate(event_ts = event_ts / 1000)

        # align arduino data to computer clock times ~~~~~~~~~~~~~~~~~~~~~
        # prepare arduino behavioral dataframe
        df_beh_events_and_tick <-
          events_arduino %>%
          mutate(event_ts = event_ts) %>%
          rename(ts = event_ts)

        # save ticks from behavior file
        ticks <- df_beh_events_and_tick %>%
          filter(event_id_char == "tick_clock_out") %>%
          group_by(session_name, blockname) %>%
          mutate(tick_number = row_number())

        # save ticks from npm
        df_rec_system_tick <- session_clock %>%
          select(session_name, blockname, ts = ts_comp) %>%
          group_by(blockname) %>%
          filter(ts > min(ts)) %>% # remove first tick from session satart
          ungroup() %>%
          mutate(event_id_char = "tick_clock_in")

        # align each block individually
        for (arduino_blockname in key_clock$blockname) {
          # perform alignment for each blockname
          events_arduino_aligned_loop <-
            convert_beh_arduino_to_fp_time(
              df_rec_system_tick %>% filter(blockname == arduino_blockname),
              df_beh_events_and_tick %>% filter(blockname == arduino_blockname)
            )

          # save session_name and blockname
          events_arduino_aligned_loop <- events_arduino_aligned_loop %>%
            mutate(session_name = session_process, blockname = arduino_blockname) %>%
            select(session_name, blockname, everything())

          if (arduino_blockname == key_clock$blockname[1]) {
            events_arduino_aligned <- events_arduino_aligned_loop
          } else {
            events_arduino_aligned <- events_arduino_aligned_loop %>%
              bind_rows(events_arduino_aligned, .)
          }
        }


        events_arduino_aligned <- events_arduino_aligned %>%
          rename(event_ts = ts_rec_system, event_ts_arduino = ts) %>%
          mutate(input_id = "arduino_aligned")

        # combine with events
        events <- events %>%
          bind_rows(events_arduino_aligned)
      }

      # save event number for all events
      events <- events %>%
        arrange(session_name, blockname, input_id, event_ts) %>%
        group_by(session_name, blockname, input_id, event_id_char) %>%
        mutate(event_number = row_number())



      # filter NAs from events
      events <- events %>%
        filter(!is.na(event_ts))

      #  ~ filter all events ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      key_events_of_interest <- log_npm %>%
        filter(session_name == session_process) %>%
        select(blockname, events_of_interest) %>%
        unique() %>%
        filter(!is.na(events_of_interest))


      # stop program and display error if there are multiple unique sets of events_of_interst for any blocks
      key_events_of_interest_counts <- key_events_of_interest %>%
        group_by(blockname) %>%
        summarise(events_of_interst_count = n(), .groups = "drop")

      if (sum(key_events_of_interest_counts$events_of_interst_count > 1) > 0) {
        print(
          "!error: more than one set of unique events_of_interst defined for following blocks"
        )
        print(key_events_of_interest_counts %>% filter(events_of_interst_count > 1))
        return(0)
      }


      # create empty tibble for events_filtered
      events_filtered <- tibble(
        session_name = character(),
        blockname = character(),
        input_id = character(),
        event_id_char = character(),
        event_ts = double(),
        event_ts_offset = double(),
        event_ts_arduino = double(),
        event_number = integer()
      )

      for (blockname_process in key_events_of_interest$blockname) {
        events_of_interest <-
          key_events_of_interest$events_of_interest[key_events_of_interest$blockname == blockname_process] %>%
          str_split(";")

        events_of_interest <- events_of_interest[[1]] %>% str_trim()

        events_filtered_loop <- events %>%
          ungroup() %>%
          filter(blockname == blockname_process) %>%
          fp_filter_events(events_of_interest)

        events_filtered <- events_filtered %>%
          bind_rows(events_filtered_loop)
      }

      print("")
      print("- processing peri-event time histograms")

      # print warning if there are no events for peri-event time histogram
      if (nrow(events_filtered) == 0) {
        print(str_c(" * no perievent time histogram: no filtered time stamps"))
      } else {
        for (blockname_process in blocknames) {
          events_filtered_loop <- events_filtered %>%
            filter(blockname == blockname_process)

          streams_data_smoothed_fitted_loop <-
            streams_data_smoothed_fitted %>%
            filter(blockname == blockname_process) %>%
            ungroup()


          if (nrow(events_filtered_loop) == 0) {
            print(str_c(" - ", blockname_process, ": no filetered time stamps"))
          } else {
            print(str_c(" - ", blockname_process, ": computing peth for ", dim(events_filtered_loop)[1], " events"))

            # re-calculate event_number
            events_filtered_loop <- events_filtered_loop %>%
              arrange(event_ts) %>%
              group_by(event_id_char) %>%
              mutate(event_number = row_number())

            log_npm_loop <- log_npm %>%
              filter(blockname == blockname_process) %>%
              select(blockname, peth_pre, peth_post, down_sampled_freq) %>%
              unique()

            # display error if peth_pre, peth_post, and/or down_sampled_freq are inconsistent within block
            if (nrow(log_npm_loop) > 1) {
              print("!error: more than one set of unique set of peth parameters in block")
              print(log_npm_loop)
              return(0)
            }

            # # # compute peth
            streams_data_smoothed_fitted_peth_downsampled_loop <-
              fp_peri_event_time_histogram(
                streams_data_smoothed_fitted_loop %>%
                  select(-session_name, -blockname_temp),
                events_filtered_loop %>%
                  select(
                    -session_name,
                    -event_ts_offset,
                    -event_ts_arduino,
                    -input_id
                  ),
                log_npm_loop$peth_pre,
                log_npm_loop$peth_post,
                log_npm_loop$down_sampled_freq,
                log_npm_loop$down_sampled_freq
              ) %>%
              ungroup()

            streams_data_smoothed_fitted_peth_downsampled_loop <-
              streams_data_smoothed_fitted_peth_downsampled_loop %>%
              mutate(session_name = session_process) %>%
              select(session_name, everything())

            if (blockname_process == blocknames[1]) {
              streams_data_smoothed_fitted_peth_downsampled <-
                streams_data_smoothed_fitted_peth_downsampled_loop
            } else {
              streams_data_smoothed_fitted_peth_downsampled <-
                streams_data_smoothed_fitted_peth_downsampled_loop %>%
                bind_rows(
                  streams_data_smoothed_fitted_peth_downsampled,
                  .
                )
            }
          }
        }
      }


      # save processed data
      print("")
      print(str_c(" - saving files to directory ", dir_npm_processed))

      for (blockname_process in blocknames) {
        print(str_c("  * blockname ", blockname_process))

        # save streams session
        print(str_c("   ~ ", blockname_process, "_streams_session.feather"))

        streams_data_smoothed_fitted %>%
          filter(blockname == blockname_process) %>%
          arrow::write_feather(
            str_c(
              dir_npm_processed,
              "/",
              blockname_process,
              "_streams_session.feather"
            )
          )

        if (nrow(events_filtered) > 0) {
          events_filtered_loop <- events_filtered %>%
            filter(blockname == blockname_process)

          if (nrow(events_filtered_loop) > 0) {
            # save events peth
            print(str_c("   ~ ", blockname_process, "_events_peth.feather"))

            events_filtered_loop %>%
              arrow::write_feather(
                str_c(
                  dir_npm_processed,
                  "/",
                  blockname_process,
                  "_events_peth.feather"
                )
              )

            # save streams peth
            print(str_c("   ~ ", blockname_process, "_streams_peth.feather"))

            streams_data_smoothed_fitted_peth_downsampled %>%
              filter(blockname == blockname_process) %>%
              arrow::write_feather(
                str_c(
                  dir_npm_processed,
                  "/",
                  blockname_process,
                  "_streams_peth.feather"
                )
              )

            # save all events
            print(str_c("   ~ ", blockname_process, "_events.feather"))

            events %>%
              filter(blockname == blockname_process) %>%
              arrow::write_feather(str_c(
                dir_npm_processed,
                "/",
                blockname_process,
                "_events.feather"
              ))
          }
        }
      }
    }

    print("")
    print("fp_preprocess() complete")
  }

remove_opto_stim_frames_npm <- function(df_fp, theshold_sd) {
  # filters out frames where red channel shows substantially higher signal compared to baseline
  # should only be applied to datasets with red imaging + red opto stimulation and NO RED INDICATOR

  opto_stim_frames <- df_fp %>%
    mutate(t = timestamp - min(timestamp)) %>%
    gather("region", "brightness", starts_with("region")) %>%
    filter(substr(region, nchar(region), nchar(region)) == "r") %>% # filter to red channel only
    group_by(frame_counter, timestamp, led_state) %>%
    summarise(brightness_mean = mean(brightness), .groups = "drop") %>%
    filter(brightness_mean > mean(brightness_mean) + theshold_sd * sd(brightness_mean)) %>%
    pull(frame_counter)

  df_fp_filtered <- df_fp %>%
    filter(!frame_counter %in% opto_stim_frames)

  return(df_fp_filtered)
}

# In Python:
# dir_extracted = r'.\examples\tdt\extracted'
# dir_processed = r'.\examples\tdt\processed'
# log_fp = r'examples\tdt\log_data_fp_tdt.csv'
# log_fp = pd.read_csv(log_fp)
# fp_preprocess(dir_extracted, dir_processed, log_fp)

# In R:
dir_extracted <- "examples/tdt/extracted"
dir_processed <- "examples/tdt/processed"
log_fp <- "examples/tdt/log_data_fp_tdt_single_session.csv"
log_fp <- read_csv(log_fp)
fp_preprocess(dir_extracted, dir_processed, log_fp, overwrite=1)