# Title     : f_general_functions.R
# Objective : collect functions for functions that are used across multiple applications
# Created by: Adam Gordon-Fennell (agg2248@uw.edu)
#             Garret Stuber Lab, University of Washington
# Created on: 1/31/2022

# general --------------------------------------------------------------------------------------------------------------
# dir_diff (used in tidy_med, fp_preprocessing)
dir_diff <- function(dir_input, dir_output, list_suffix, print_message = 0, ignore_suffixes = NA){
  # return the files located in dir_input that are not already located in dir_output
  #
  # list_suffix provides strings that contain file name suffix + extension to be removed off of file name

  list_dir_input <- list.files(dir_input)
  list_dir_output <- list.files(dir_output)

  # filter out files in dir_output that contain suffix to ignore
  if(sum(!is.na(ignore_suffixes)) >= 1){
    for(ignore_suffix in ignore_suffixes){
      list_dir_output <- list_dir_output[!str_detect(list_dir_output, ignore_suffix)]
    }
  }

  for(suffix in list_suffix){
    list_dir_input <- list_dir_input %>%
      str_remove(suffix)

    list_dir_output <- list_dir_output %>%
      str_remove(suffix)
  }

  list_dir_input <- list_dir_input %>% unique()
  list_dir_output <- list_dir_output %>% unique()

  # if there are any files in the output directory that are not in the input directory, give a warning
  if(length(setdiff(list_dir_output, list_dir_input))){
    if(print_message){
      print("WARNING: THE FOLLOWING FILES ARE LOCATED IN THE dir_output THAT ARE NOT FOUND IN dir_input")
    }

    for(fn_missing in setdiff(list_dir_output, list_dir_input)){
      if(print_message){
        print(fn_missing)
      }
    }
    print("")
  }

  # create a list of files contained in dir_input but not dir_output
  fns <- setdiff(list_dir_input, list_dir_output)

  if(print_message){
    print(str_c("number of files not found in dir_output: ", length(fns)))
    print("")
  }

  return(fns)
}

dir_diff_session <- function(dir_input, dir_output, list_suffix, print_message = 0, log_data){
  # return the files located in dir_input that are not already located in dir_output
  #
  # list_suffix provides strings that contain file name suffix + extension to be removed off of file name

  list_dir_input <- list.files(dir_input)
  list_dir_output <- list.files(dir_output)

  for(suffix in list_suffix){
    list_dir_input <- list_dir_input %>%
      str_remove(suffix)

    list_dir_output <- list_dir_output %>%
      str_remove(suffix)
  }

  list_dir_input <- list_dir_input %>% unique()
  list_dir_output <- list_dir_output %>% unique()


  list_blockname <- log_data %>%
    filter(session_name %in% list_dir_input) %>%
    pull(blockname) %>%
    unique()

  blocknames <- setdiff(list_blockname, list_dir_output)

  session_names <- log_data %>%
    filter(blockname %in% blocknames) %>%
    pull(session_name)

  if(print_message){
    print(str_c("number of files not found in dir_output: ", length(fns)))
    print("")
  }

  return(session_names)
}

# generic combine / append function
combine_files <- function(dir_input_tidy, fns, output_fn, suffix, append_combined){
  # combines multiple files into a single combined file. If the combined file already exists, new files will be appended
  # to the existing combined file
  #
  # dir_input_tidy: directory of input tidy data
  # fns: list of files within dir_input_tidy that we will combine
  # suffix: suffix of input files and output file (e.g. "_events.csv")
  # append_combined: 1: append to existing combined dataframe, 0: genreate from scratch and overwrite existing dataframe

  skip <- 0
   if(append_combined){
      if(file.exists(output_fn)){
        if(str_detect(suffix, 'csv')){fns_combined <- read_csv(output_fn, col_types = cols())}
        if(str_detect(suffix, 'feather')){fns_combined <- read_feather(output_fn)}

        fns_combined <- fns_combined %>% pull(med_file_name) %>% unique()

        fns_diff <- setdiff(fns %>% str_remove(suffix), fns_combined)

        if(length(fns_diff) > 0){
          fns_diff <- fns_diff %>%
            str_c(suffix)
        }

        if(length(fns_diff) == 0){
          print(str_c("all raw files already found in combined dataframe ", output_fn))
          skip <- 1
        } else {
          fns <- fns_diff
        }
      } else {
        fns_diff <- vector('character')
      }
    } else {
      print("append overwrite... creating new output dataframe")
      fns_diff <- vector('character')
    }

    if(length(fns) == 0){
      print(str_c("WARNING NO FILES FOUND WITH fn_string = ", fn_string, " & suffix = ", suffix))
      next
    }

  if(skip != 1){
      print(str_c("combining ", length(fns), " files with suffix ", suffix))

      first_file <- 1

      for(fn in fns){
        if(str_detect(fn, 'csv')){df_loop <- read_csv(str_c(dir_input_tidy, '/', fn), col_types = cols())}
        if(str_detect(fn, 'feather')){df_loop <- read_feather(str_c(dir_input_tidy, '/', fn))}
        if(!str_detect(fn, 'csv') & !str_detect(fn, 'feather')){
          print(str_c("WARNING INCOMPATABLIE FILE TYPE FOR FILE ", fn))
          next
        }

        if(first_file == 1){
          df_out <- df_loop
          first_file <- 0
        } else {
          df_out <- df_loop %>% bind_rows(df_out, .)
        }
      }

      print(str_c("writing: ", output_fn))

      if(str_detect(suffix, 'csv')){
        if(append_combined & length(fns_diff) > 0){
          print(output_fn)
          print(read_csv(output_fn))
          print(df_out)
          df_out %>% bind_rows(read_csv(output_fn, col_types = cols()),.) %>% write_csv(output_fn)
        } else {
          df_out %>% write_csv(output_fn)
        }
      }

      if(str_detect(suffix, 'feather')){
        if(append_combined & length(fns_diff) > 0){
          df_out %>% bind_rows(read_feather(output_fn, col_types = cols()),.) %>% write_feather(output_fn)
        } else {
          df_out %>% write_feather(output_fn)
        }
      }

      print("")
    }
  }

coerce_character <- function(t1, t2){
  # compare the datatypes of 2 dataframes and then convert numeric data to character if there is a disagreement
  #
  t1_dtype <- t1 %>%
    head() %>%
    collect() %>%
    lapply(class) %>%
    unlist() %>%
    as_tibble() %>%
    filter(!value %in% c('difftime', 'POSIXt')) %>% # filter out double data type for dates
    rename(t1_dtype = value) %>%
    mutate(variable = colnames(t1))


  t2_dtype <- t2 %>%
    head() %>%
    collect() %>%
    lapply(class) %>%
    unlist() %>%
    as_tibble() %>%
    filter(!value %in% c('difftime', 'POSIXt')) %>% # filter out double data type for dates
    rename(t2_dtype = value) %>%
    mutate(variable = colnames(t2))

  diff_t1_t2 <-
    suppressMessages(left_join(t1_dtype, t2_dtype)) %>%
    filter(t1_dtype != t2_dtype)

  if(nrow(diff_t1_t2 > 0)){
    for(n_var in 1:nrow(diff_t1_t2)){
      var_name <- diff_t1_t2$variable[n_var]
      if(diff_t1_t2$t1_dtype[n_var] == 'character'){
        t2 <- t2 %>% mutate_at(vars(var_name), as.character)
      }

      if(diff_t1_t2$t2_dtype[n_var] == 'character'){
        t1 <- t1 %>% mutate_at(vars(var_name), as.character)
      }
    }
  }

  return(list(t1, t2))
}



combined_import <- function(import_directory, filter_strings = NA, prefixes = NA, suffixes, verbose = 0){
  # combine a set of csv or feather files into a single dataframe
  #
  # inputs:
  # - import_directory: directory of files
  # - filter_strings (OPTIONAL): vector of strings to filter files in import_directory with
  # - prefix (OPTIONAL): vector of strings to filter files in import_directory with
  # - suffixes: vector of strings of suffixes you wish to combine
  #
  # output:
  # - list of dataframes that combine files with specified suffixes

    if(missing(verbose)) {
        verbose <- 0
    }

  if(verbose){
    print('combined import...')
    print('')
  }

  if(length(filter_strings) == 0){
    print('combined_import cancled: filter_strings has length of 0')
    return(0)
  }

  num_suffix <- 1

  for(suffix in suffixes){
    fns <- list.files(import_directory)

    if(verbose){
      print(str_c('number of files in ', import_directory, ': ', length(fns)))
      print('')
    }

    # filter fns to files with  strings included in filter_strings
    if(sum(!is.na(filter_strings)) > 0){
      for(filter_string in filter_strings){
        fns_loop <- fns[str_detect(fns, filter_string)]

        if(filter_string == filter_strings[1]){
          fns_combined <- fns_loop
        } else{
          fns_combined <- c(fns_combined, fns_loop)
        }
      }
    fns <- fns_combined
    }

    # filter fns to files with  strings included in prefixes
    if(sum(!is.na(prefixes)) > 0){
      for(prefix in prefixes){
        fns_loop <- fns[str_detect(fns, prefix)]

        if(prefix == prefixes[1]){
          fns_combined <- fns_loop
        } else{
          fns_combined <- c(fns_combined, fns_loop)
        }
      }
    fns <- fns_combined
    }

    # filter fns to files with suffix
    fns <- fns[str_detect(fns, suffix)]

    if(length(fns) == 0){
      print(str_c('No file names matching criteria found in directory ', import_directory))
      print('')
      print('check suffixes and extension')
      return(0)
    }

    if(verbose){
      print(str_c('combining ', length(fns), ' files'))
      print('')
    }

    # read in and combine files in fitler_list
    if(str_detect(suffix, 'feather')){
      for(fn in fns){
        if(verbose){print(fn)}

        data_loop <- read_feather(str_c(import_directory, '/', fn)) %>%
          mutate(fn = fn) %>%
          select(fn, everything())

        if(fn == fns[1]){
          data_combined <- data_loop
        } else {
          data_combined_list <- coerce_character(data_combined, data_loop)

          data_combined <- data_combined_list[[1]]
          data_loop     <- data_combined_list[[2]]

          data_combined <- data_loop %>% bind_rows(data_combined,.)
        }
      }
    }

    if(str_detect(suffix, 'csv')){
      for(fn in fns){
        if(verbose){print(fn)}

        data_loop <- read_csv(str_c(import_directory, '/', fn), col_types = cols()) %>%
          mutate(fn = fn) %>%
          select(fn, everything())

        if(fn == fns[1]){
          data_combined <- data_loop
        } else {
          data_combined_list <- coerce_character(data_combined, data_loop)

          data_combined <- data_combined_list[[1]]
          data_loop     <- data_combined_list[[2]]

          data_combined <- data_loop %>% bind_rows(data_combined,.)
        }
      }
    }

    if(num_suffix == 1){
      output_list <- list(data_combined)
      } else {
        output_list[[num_suffix]] <- data_combined
        }

    num_suffix <- num_suffix + 1
  }

  return(output_list)
}


# combine based on string in file name
combine_dfs_fn_string <- function(dir_input_tidy, dir_output_combined, list_suffix, fn_string, append_combined){
  # combine multiple tidy datasets into a single file
  #
  # inputs:
  # - dir_input_tidy (string): system directory of tidy datasets you wish to combine
  # - dir_output_combined (string): system directory of combined tidy datasets will be saved to
  # - list_suffix (vector of strings): file name suffixes contained in dir_input_tidy you wish to combine
  #   must include file extensions .csv or .feather
  #   NOTE: the input suffix will be used to set the output data type (e.g. .csv in .csv out)
  # - fn_string (string): string contained in file names within dir_input_tidy you wish to combine into a single file
  # - append_combined (boolean): 1: append to existing combined dataframe, 0: genreate from scratch and overwrite
  #   existing dataframe
  #
  # outputs:
  #  - combined tidy datasets (1 / member of list_suffix)

  list_dir_input_tidy_files <- list.files(dir_input_tidy)

  list_dir_input_tidy_files <- list_dir_input_tidy_files[list_dir_input_tidy_files %>% str_detect(fn_string)]

  for(suffix in list_suffix){
    fns <- list_dir_input_tidy_files[list_dir_input_tidy_files %>% str_detect(suffix)]

    output_fn <- str_c(dir_output_combined, '/data_med_combined_', fn_string, suffix)

    combine_files(dir_input_tidy, fns, output_fn, suffix, append_combined)
  }
}

combine_dfs_log_variable <- function(dir_input_tidy, dir_output_combined, list_suffix, log_data, variables, values, append_combined){
  # combine multiple tidy datasets into a single file
  #
  # inputs:
  # - dir_input_tidy (string): system directory of tidy datasets you wish to combine
  # - dir_output_combined (string): system directory of combined tidy datasets will be saved to
  # - list_suffix (vector of strings): file name suffixes contained in dir_input_tidy you wish to combine
  #   must include file extensions .csv or .feather
  #   NOTE: the input suffix will be used to set the output data type (e.g. .csv in .csv out)
  # - log_data (dataframe): dataframe containing a variable with file names, and variables used to define data to combine
  # - variables (vector of strings): variable names contained in dataframe you wish to combine
  # - values (vector of strings): values of variables that you wish to filter to
  # - append_combined (boolean): 1: append to existing combined dataframe, 0: genreate from scratch and overwrite
  #   existing dataframe
  #
  # outputs:
  #  - combined tidy datasets (1 / member of list_suffix)

  # filter log_data that have variable(n) that matches value(n)
  for(n_filter in 1:length(variables)){
    log_data <- log_data %>%
      filter((!!sym(variables[n_filter])) == values[n_filter])
  }

  list_dir_input_tidy_files <- list.files(dir_tidy)

  for(suffix in list_suffix){
    list_dir_input_tidy_files <- list_dir_input_tidy_files %>%
      str_remove(suffix) %>%
      unique()

  }

  fn_stems <- setdiff(list_dir_input_tidy_files, log %>% pull(med_file_name))

  for(suffix in list_suffix){
    # generate file names
    fns <- fn_stems %>% str_c(suffix)

    # create output_fn
    output_fn <- str_c(dir_output_combined, '/data_med_combined')

    for(n_filter in 1:length(variables)){
      output_fn <- output_fn %>% str_c('__', variables[n_filter], '_', values[n_filter])
    }

    output_fn <- output_fn %>% str_c(suffix)

    combine_files(dir_input_tidy, fns, output_fn, suffix, append_combined)
  }
}


convert_day_character_to_numeric <- function(df, ...){
  # converts character strings for days to integers relative to first day (e.g. 'd01' -> 1)
  df %>%
    mutate(day = day %>% str_remove('d') %>% as.integer()) %>%
    group_by_(...) %>%
    mutate(day = day - min(day) + 1) %>%
    ungroup()
}

get_day_numeric <- function(df, ...){
  # converts character strings for days to integers relative to first day (e.g. 'd01' -> 1)

  if('day' %in% colnames(df)){
    print('error: variable "day" already in dataframe')
  } else {
  df %>%
    select(subject, date) %>%
    unique() %>%
    arrange(subject, date) %>%
    group_by(subject) %>%
    mutate(day = row_number()) %>%
    left_join(df,.)
  }
}

unique_values <- function(df,...){
  df %>%
    ungroup() %>%
    select(...) %>%
    unique()
}

format_dir <- function(x){
  if(str_sub(x, nchar(x), nchar(x)) != '/'){
    x <- x %>% str_c('/')
  }
  return(x)
}


get_recent_para_session <- function(dir_para_session) {
  # return most recent log_data and key_pump from para_session_yyyy_mm_dd.xlsx file based on date

  dir_para_session <- format_dir(dir_para_session)
  fn_parasessions <- list.files(dir_para_session)

  fn_parasessions_date <- fn_parasessions %>%
    str_remove('para_session_') %>%
    str_remove('.xlsx') %>%
    ymd()

  fn_parasession_most_recent <- fn_parasessions[fn_parasessions_date == max(fn_parasessions_date)]

  log_data <- read_xlsx(str_c(dir_para_session, fn_parasession_most_recent), na = c('','NA')) %>%
    mutate(date = ymd(date))

  key_pump <- read_xlsx(str_c(dir_para_session, fn_parasession_most_recent), sheet = 'pump_index', na = c('','NA')) %>%
    mutate(date = ymd(date))

  print(str_c('returning most recent para_session: ', fn_parasession_most_recent))
  return(list(log_data,key_pump))

}

# Get events -----------------------------------------------------------------------------------------------------------
get_event_bouts <- function(events, event_id_char_of_interest,  filt_tm, filt_n, filt_dir, id_char){
  # extract either onset or offsets of bouts of a chosen event based on the number of events prior/post (filt_n) within
  # a timeframe prior/post (filt_dir)
  #
  # inputs:
  #  - events (df): dataset from a single session that includes event_id_char (the id of the event) and event_ts (the
  #    time of the event
  #  - event_id_char_of_ineterest (char vector): vector of strings that contains the events that you want to use
  #  - filt_tm (double vector): time for window prior and post event that will be used
  #      - Units of this variable must match units for event_ts in events
  #  - filt_n (integer vector): number of events in the window prior and post that will be used
  #  - filt_dir(character vector): logical used to apply to time window ('>', '<', or '==')
  #  - id_char (string): string to replace event_id_char in output dataframe
  #
  # directions:
  #  - filt_tm, filt_n, and filt_dir each contain 2 values that apply to the window prior to and post each event
  #  - events will be filtered down to events that match the specified conditions
  #
  # example:
  #  - licking bout onsets
  #  - condition: no licks in 1s prior to first lick, at least 1 lick in 1s post lick
  #  input: events_get_event_bouts(events, c('lick'), c(1,1), c(0, 0), c('==', '>'), 'lick_onset')

  # filter to data that match the event of interest
  events_bout_onset <- events %>%
    filter(event_id_char %in% event_id_char_of_interest) %>%
    mutate(n_event_prior = 0,
           n_event_post  = 0) %>%
    arrange(event_ts)

  event_tss <- events_bout_onset %>%
    pull(event_ts)

  # for each time stamp, retrieve the number of events preceding and following
  for (event in 1:nrow(events_bout_onset)){
    events_bout_onset$n_event_prior[event] <-
      sum(
        event_tss > events_bout_onset$event_ts[event] - filt_tm[1] &
        event_tss < events_bout_onset$event_ts[event]
        )

    events_bout_onset$n_event_post[event] <-
      sum(
        event_tss > events_bout_onset$event_ts[event] &
        event_tss < events_bout_onset$event_ts[event] + filt_tm[2]
        )
  }

  # filter based on events prior
  if(filt_dir[1] == '<'){
    events_bout_onset <- events_bout_onset %>%
    filter(n_event_prior < filt_n[1])
  }

  if(filt_dir[1] == '>'){
    events_bout_onset <- events_bout_onset %>%
    filter(n_event_prior > filt_n[1])
  }

  if(filt_dir[1] == '=='){
    events_bout_onset <- events_bout_onset %>%
    filter(n_event_prior == filt_n[1])
  }

  # filter based on events post
  if(filt_dir[2] == '<'){
    events_bout_onset <- events_bout_onset %>%
    filter(n_event_post < filt_n[2])
  }

  if(filt_dir[2] == '>'){
    events_bout_onset <- events_bout_onset %>%
    filter(n_event_post > filt_n[2])
  }

  if(filt_dir[2] == '=='){
    events_bout_onset <- events_bout_onset %>%
    filter(n_event_post == filt_n[2])
  }

  # rename event_id_char so it can be combined with input dataset
  events_bout_onset <- events_bout_onset %>%
    mutate(event_id_char = id_char) %>%
    select(-n_event_prior, -n_event_post)

  return(events_bout_onset)
}

# time series functions ------------------------------------------------------------------------------------------------

# create perievent dataframe from series and event dataframes
return_peri_event_series <- function(df_series, df_events, var_series_ts, var_event_ts, var_event_num, var_grouping, time_pre, time_post, filt_duplicate){
  # return a perievent time series with relative timing
  #
  # inputs:
  # - df_series: dataframe that contains one or more time series in tidy format
  # - df_events: dataframe that contains the timing of events (should match sampling rate of df_series)
  # - var_series_ts: string of the variable name for time within df_series
  # - var_event_ts:  string of the variable name for time within df_events (CANNOT NOT MATCH var_series_ts!)
  # - var_event_num: string of the variable name for the event number within df_events
  # - var_grouping:  string of the variable names(s) within df_series and df_events that should be used to filter df_series
  #     e.g. multiple time series are contained within df_series and df_events (fiber A / B, day 1 / 2)
  # - time_pre: time prior to event that should be contained in output dataframe (units match inputs)
  # - time_post: time following event that should be contained in output dataframe (units match inputs)
  # - filt_duplicate: boolean- 1: filter out subsequent events that occured within time_post

  # output:
  # - peri_event_series
  #   will contain all columns found within df_series filtered based on var_event_ts and var_series_ts
  #   will have added the var_event_ts and var_event_num contained in df_events
  #   will also contain rel_time (df_series$var_series_ts - df_events$var_event_ts)
  #
  #   note: rel_time may need to be rounded

  event_time_loop_previous <- NA

  for(n_event in 1:nrow(df_events)){
    loop_event <- df_events[n_event,] # save row from events dataframe for loop

    loop_series <- df_series # setup series dataframe for loop


    event_time_loop <- loop_event %>% pull(var_event_ts)


    if(filt_duplicate){
      if(!is.na(event_time_loop_previous)){
        if(event_time_loop - event_time_loop_previous < time_post & event_time_loop - event_time_loop_previous > 0 ){
          next # go to next loop
        }
      }
    }

    # filter the series data frame using the loop_event based on variables defined in var_grouping
    for(var in var_grouping){
      loop_series <- loop_series %>%
        filter((!!as.symbol(var)) == loop_event %>% pull(var))
    }

    # filter series based on event time in loop_event
    loop_series <- loop_series %>%
      filter((!!as.symbol(var_series_ts)) > loop_event %>% pull(var_event_ts) - time_pre,
             (!!as.symbol(var_series_ts)) < loop_event %>% pull(var_event_ts) + time_post)

    # save loop_event info in series
    loop_series <- loop_series %>%
      mutate(!!var_event_ts  := loop_event %>% pull(var_event_ts),
             !!var_event_num := loop_event %>% pull(var_event_num)
             )

    # create peri event time
    loop_series <- loop_series %>%
      mutate(rel_time := !!as.name(var_series_ts) - !!as.name(var_event_ts))

    if(n_event == 1){
      peri_event_series <- loop_series
    } else {
      peri_event_series <- loop_series %>% bind_rows(peri_event_series,.)
    }

    event_time_loop_previous <- event_time_loop
  }

  return(peri_event_series)
}

# compute binned average of perievent data
return_event_binned_counts <- function(df_events, var_relative_ts, time_start, time_end, bin_width, ...){
  # return the binned counts for events within a dataframe
  #
  # inputs:
  # - df_events: dataframe that contains the relative time of events
  # - var_relative_ts: string of the variable name for relative time within df_peri_event
  # - time_start: time start (units match var_relative_ts in df_peri_event)
  # - time_end: time followig event that should be included in output dataframe (units match var_relative_ts in df_peri_event)
  # - bin_width: the width of time bins for computing binned counts (units match var_relative_ts in df_peri_event)
  # - ...: grouping variables
  #
  # output:
  # - peri_event_binned: dataframe with following varibles
  #   - grouping variables as defined by var_grouping
  #   - event_bin: left side of time bin
  #   - event_count: number of observations that occured during the corresponding bin
  var_relative_ts <- enquo(var_relative_ts)

  bins <- seq(time_start, time_end, bin_width)

  df_event_binned <- df_events %>%
    mutate(event_bin = cut(!!var_relative_ts, bins, label = bins[1:length(bins) - 1])) %>%
    mutate(event_bin = event_bin %>% as.character() %>% as.double()) %>%
    mutate(event_bin = event_bin %>% round(4) %>% as.character()) %>% # need to convert back to character to avoid slight numerical differences
    group_by(..., event_bin) %>%
    summarise(event_count = n()) %>%
    ungroup() %>%
    complete(..., event_bin = as.character(bins[1:length(bins) - 1]), fill = list(event_count = 0)) %>%
    mutate(event_bin = event_bin %>% as.double())

  return(df_event_binned)
}

get_peaks_threshold <- function(streams, var_time, var_dv, fs, filt_window, var_stream_id, peak_threshold){
  # return peak times from streams
  #
  # inputs
  #  streams (dataframe): time series data
  #  var_time(string): variable name of time within streams
  #  var_dv (string): variable name of signal within streams
  #  fs (double): sampling rate (Hz)
  #  filt_window(double): time window to filter double catches (set to 0 to return all peak times)
  #  var_stream_id(vector of strings): grouping variable names (e.g. c('blockname', 'channel_id'), or c('etl_plane', 'cell_number'))
  #  peak_threshold(double): threshold for peak detection (units match units of var_dv)

  streams <- streams  %>%
    rename(time = !!as.name(var_time))

  # return variable that denotes if sample is above peak_threshold
  streams <- streams %>%
    arrange_(c(var_stream_id, 'time')) %>%
    group_by(.dots = lapply(c(var_stream_id), as.symbol)) %>%
    mutate(above_threshold_window = ifelse(!!as.name(var_dv) >= peak_threshold,  1, 0)) %>%
    ungroup()

  # return peak onset times
  temp_peaks <- streams %>%
    group_by(.dots = lapply(c(var_stream_id), as.symbol)) %>%
    filter(above_threshold_window == 1) %>% # filter to samples above threshold
    mutate(time_step = time - lag(time)) %>% # compute time step
    mutate(onset_above_threshold_window  = ifelse(time_step > 1/fs + 1/(fs*10) | is.na(time_step), 1, NA)) %>% # filter to samples are not sequential
    filter(!is.na(onset_above_threshold_window)) %>%
    select(-time_step) %>%
    mutate(peak_num = row_number()) %>%
    ungroup()

  temp_peaks <- suppressMessages(left_join(streams,temp_peaks))


  # return peak points within periods above peak_threshold
  temp_peaks <- temp_peaks %>%
    group_by(.dots = lapply(c(var_stream_id), as.symbol)) %>%
    fill(peak_num) %>%
    group_by(.dots = lapply(c(var_stream_id, 'peak_num'), as.symbol)) %>%
    mutate(peak = ifelse(!!as.name(var_dv) == max(!!as.name(var_dv)), 1, 0)) %>% #
    ungroup()


  # filter peak times
  temp_peaks_filt <- temp_peaks %>%
    filter(peak == 1) %>%
    rename(event_ts = time) %>%
    mutate(event_id_char = 'peak') %>%
    group_by(.dots = lapply(c(var_stream_id), as.symbol)) %>%
    do(get_event_bouts(.,
                       event_id_char_of_interest = c('peak'),
                       filt_tm = c(filt_window, filt_window),
                       filt_n = c(0, 9999),
                       filt_dir = c("==", "<"),
                       id_char = 'filtered_peak'
                       )
       )%>%
      rename(time = event_ts) %>%
      select(var_stream_id, time, event_id_char) %>%
      left_join(temp_peaks,., by = c(var_stream_id, 'time')) %>%
    filter(event_id_char == 'filtered_peak') %>%
    group_by(.dots = lapply(c(var_stream_id), as.symbol)) %>%
    mutate(peak_num = row_number()) %>%
    rename(event_time = time) %>%
    select(-above_threshold_window, -onset_above_threshold_window, -peak, -event_id_char)

  return(temp_peaks_filt)
}

# misc pre processing functions
return_trial_ids <- function(dir_extracted, manual_fns, file_format_output){
  # returns trial_ids for trial based data collected with arduino system
  #  trial based vbariables should begin with prefix "trial_"

  data_list <- combined_import(import_directory = dir_extracted,
                               filter_strings = manual_fns,
                               prefixes = NA,
                               suffixes = c(str_c('_param_dynamic.', file_format_output))
                               )

  param_dynamic <- data_list[[1]]

  param_dynamic <- param_dynamic %>%
    filter(param_dynamic %>% str_detect('trial_')) %>%
    arrange(blockname, param_ts, param_dynamic) %>%
    select(blockname, param_dynamic, param_value) %>%

    group_by(blockname, param_dynamic) %>%
    mutate(trial_num = row_number()) %>%
    spread(param_dynamic, param_value) %>%
    group_by(blockname) %>%
    filter(trial_num < max(trial_num))

  return(param_dynamic)
}

# finalizing analysis --------------------------------------------------------------------------------------------------
write_data_local <- function(dir_raw, dir_local, blocknames, extension){
  blocknames <- str_c(blocknames, '.', extension)

  if(dir.exists(dir_raw)){
    for(blockname in blocknames){
      file_path <- str_c(format_dir(dir_raw), blockname)

      if(file.exists(file_path)){

        read_csv(file_path) %>%
          write_csv(str_c(format_dir(dir_local), blockname))

      } else {
        print(str_c('file: ', file_path, ' does not exist'))
      }
    }
  } else {
    print(str_c('directory: ', dir_raw, ' does not exist'))
  }

}
