# Title     : f_med_associates.R
# Objective : collect functions for tidying, combining, and analyzing med associates data
# Created by: Adam Gordon-Fennell (agg2248@uw.edu)
#             Garret Stuber Lab, University of Washington
# Created on: 1/31/2022

# tidy_med--------------------------------------------------------------------------------------------------------------

# line_getkey (used in tidy_med)
line_getkey <- function(ln){
  if(nchar(ln) < 2){
    return(FALSE)
  } else if(substr(ln,1,1) == ' '){
    return(FALSE)
  } else if(grepl(':', ln)){
    return(str_split(ln,':')[[1]][1])
  } else {
    return(FALSE)
  }
}

# line_getval (used in tidy_med)
line_getval <- function(ln){
  if(nchar(ln) < 2){
    return(ln)
  } else if(grepl(':', ln)){
    parts <- strsplit(ln, ':')[[1]]
    return(parts[2:length(parts)] %>% str_c(collapse = ":"))
  } else {
    return(ln)
  }
}

tidy_med <- function(dir_input, dir_output, key_med_variables, key_event_ids, log_data, extension, override){
  # takes raw med files and converts them into tidy csv or feather files defined in log_data
  # 
  # input: 
  #  - dir_input (string): absolute directory of raw med fiels
  #      - Sll files within the input directory will be converted.
  #  - dir_output (string): absolute directory to export csv or feather files to
  #  - key_med_variables (dataframe): key that contains the identity of variables within med file containing the
  #    following variables
  #       - key_variable_med (int): numerical identity of med key to apply to raw med data
  #       - id (string): identity of med variable
  #       - med_var (char): character identity of variable vectors within med (e.g. Q1 in med associates has a Q here)
  #       - index (int or char):
  #         - numerical index of variable within med_var (e.g. Q1 in med associates has a 1 here)
  #         - use "vector" for event time stamps and ids (e.g. one med associates vector Y for tracking events ids and
  #           a second vector Z for tracking event time stamps). Additional vectors will be stored as additional
  #           variables in the output data
  #        - example: https://docs.google.com/spreadsheets/d/1msVV9vabZQCwRRuj0HAdJWEhH-NThPrKAkOc1KjgMls
  #  - key_event_ids (dataframe): key that contains the identity of event codes within the raw med file containing the
  #    following variables
  #       - key_event_med (int): numerical identity of med key to apply to raw med data
  #       - event_id (int): nuemerical identity of events within the raw med associates file
  #       - event_id_char (string): character identiy of event_id's found within raw med file (e.g. active response
  #         recorded with event value 3 could be named "active_response" in event_id_char)
  #       - example: https://docs.google.com/spreadsheets/d/1msVV9vabZQCwRRuj0HAdJWEhH-NThPrKAkOc1KjgMls
  #  - log_data (dataframe): log of med associates data with 1 observation per behavioral session. Must contain the
  #    following variables
  #    - med_file_name (string): name of the raw med file (used to join in key_event_med)
  #    - key_variable_med (int): numerical identity of med key for variables to apply to raw med data
  #    - key_event_med (int): numerical identity of med key for events to apply to raw med data
  #  - override (boolean):
  #     0 extracts only data found in dir_input not found in dir_output
  #     1: extract all data and overwrites data in dir_output

  #
  # output:
  #   - all raw med files will have 2 corresponding output files named with the raw file name folowed by *_events or
  #     *_univariate with the chosen file extension (csv or feather)
  #
  #   - *_events.extension----------------------------------------------------------------------------------------------
  #        - this file will contain the events defined using vector within the index of the key_event_med dataframe
  #        - observation: event
  #        - variables:
  #            - n_sbj: the number of the subject contained in the raw file name
  #            - med_file_name: name of the raw med file
  #            - experiment: experiment as defined in raw med file
  #            - group: group as defined in raw med file
  #            - date: start date from raw med file
  #            - subject: subject as defined in raw med file
  #            - key_event_med: med key used to process raw data frame (defined in log_data and key_med_variables)
  #            - event_number (int): number of event within raw med file
  #            - event_id (int): numerical identity of events recorded in med associates event_id variable
  #            - event_ts (double): time stamp of event recorded in med associates event_ts variable
  #              - NOTE: the units for this variable is determined based on the med associates program
  #            - event_id_char (string): character identity of event_id as defined in key_event_ids.
  #              - NOTE: event ids that are found in raw med file but not in key_event_ids will show up as NA
  #
  #   - *_univariate.extension------------------------------------------------------------------------------------------
  #     - this file will contain the univariate data defined by a variable character and index within key_event_med
  #       - observation: session (1 row for each session contained in raw med file)
  #       - variables:
  #            - n_sbj: the number of the subject contained in the raw file name
  #            - med_file_name: name of the raw med file
  #            - experiment: experiment as defined in raw med file (if it exists in raw med file)
  #            - group: group as defined in raw med file (if it exists in raw med file)
  #            - date: start date from raw med file (if it exists in raw med file)
  #            - subject: subject as defined in raw med file (if it exists in raw med file)
  #            - one variable for each variable defined in key_event_med
  #              - column names are defined by key_med_variables id
  #              - data source is defined by key_med_variables med_var and index

  require(tidyverse)
  require(lubridate)

  if(grepl('feather', extension)){
  require(arrow)
  }

  print("extracing med associates data from")
  print(str_c("dir_input: ", dir_input))
  print("to")
  print(str_c("dir_output: ", dir_output))
  print("")

  if(override){
    fns <- list.files(dir_input)
    print("override... extracting all data in dir_input")
    print("")
  } else {
    fns <- dir_diff(dir_input, dir_output, c(str_c("_event.", extension), str_c("_univariate.", extension)))
  }

  if(length(fns) == 0){
    print("all data in dir_input is already extracted")
  }

  for (fn in fns){
    if(fn %in% (log_data %>% pull(med_file_name))){ # if med_file_name is found in log_data
      print(str_c("extracting file: ", fn))
    } else {
      print(str_c('WARNING! SKIPPING FILE NOT FOUND IN LOG_DATA: ', fn))
      next
    }

    lns <- vector()

    con = file(str_c(dir_input,  '/', fn), "r")
    n_line <- 1
    while ( TRUE ) {
      line = readLines(con, n = 1)
      if ( length(line) == 0 ) {
        break
      }
      lns[n_line] <- line
      n_line <- n_line + 1
    }

    close(con)

    # determine index for the start of each subject within med file
    index_sbj <- vector()

    for (n_line in seq(1, length(lns))){
      if (grepl('Start Date:', lns[n_line])){
        index_sbj <- append(index_sbj, n_line)
      }
    }

    data <- tibble('n_sbj' = character(),
                   'med_var' = character(),
                   'value' = character())

    for(n_sbj in seq(1,length(index_sbj))){
      current_key <- FALSE
      med_var_vec <- c()
      value_vec <- c()

      if(n_sbj < length(index_sbj)){
        loop_lns <- lns[index_sbj[n_sbj]:index_sbj[n_sbj + 1]-1]
      } else {
        loop_lns <- lns[index_sbj[n_sbj]:length(lns)]
      }

      for(ln in loop_lns){
        line_key_check <- line_getkey(ln)
        if(line_key_check != FALSE){
          line_key <- line_key_check
        }

        line_val <- line_getval(ln) %>% str_trim(side = 'left')
        line_val <- line_val %>% str_split(' ')

        if(sum(!is.na(line_val[[1]])) != 0){
          for(i_value in seq(1, length(line_val[[1]]))){
            if(line_val[[1]][i_value] != ''){
                med_var_vec <- append(med_var_vec, line_key)
                value_vec   <- append(value_vec, line_val[[1]][i_value])
            }
          }
        }
      } # end for loop_lns

      data <- data %>%
        bind_rows(tibble('med_var' = med_var_vec,
                         'value'   = value_vec) %>%
                    mutate(n_sbj = n_sbj %>% as.character()))
    }

    data <- data %>%
      mutate(med_file_name = fn)

    join_id <- data %>%
      filter(med_var %in% c('Experiment', 'Start Date', 'Subject', 'Group')) %>%
      spread(med_var, value)

    # join in id info and log_data
    data <- data %>%
      left_join(join_id, by = c("n_sbj", "med_file_name")) %>% # join session id info
      left_join(log_data %>% select(med_file_name, key_variable_med, key_event_med) %>% unique(), by = "med_file_name")  # join in med keys

    data_multivariate <- data  %>%
      left_join(key_med_variables, by = c("med_var", "key_variable_med")) %>% # join in med variables using key_variable_med
      filter(index == 'vector') %>%
      group_by(med_file_name, n_sbj, med_var) %>%
      mutate(event_num = row_number()) %>%
      ungroup()

    data_univariate <- data %>%
      anti_join(data_multivariate %>% select(n_sbj, med_var) %>% unique(), by = c("n_sbj", "med_var")) %>%
      group_by(n_sbj, med_var) %>%
      mutate(index = row_number() %>% as.character()) %>%
      left_join(key_med_variables, by = c("med_var", "key_variable_med", "index")) %>%
      filter(!is.na(id)) %>%
      ungroup() %>%
      select(n_sbj, med_file_name, key_variable_med, key_event_med, id, value) %>%
      spread(id, value) %>%
      left_join(join_id, by = c("n_sbj", "med_file_name"))

    # spread event time stamps / ids
    data_multivariate <- data_multivariate %>%
      select(-med_var) %>%
      spread(id, value)

    # join in event_id_char
    if('event_id' %in% colnames(data_multivariate)){
          data_multivariate <- data_multivariate %>%
      mutate(event_id = as.integer(event_id)) %>%
      select(-index) %>%
      left_join(key_event_ids, by = c("key_event_med", "event_id"))
    } else {
      print("event_ts not found... check to ensure that key_event_med matches across files and med_var is correct")
      next
    }



    # clean up variable names
      data_multivariate <- data_multivariate %>% clean_names()
      data_univariate <- data_univariate %>% clean_names()

    if(grepl('csv', extension)){
      data_multivariate %>%
        write_csv(str_c(dir_output, '/', fn, '_event.csv' ))

      data_univariate %>%
        write_csv(str_c(dir_output, '/', fn, '_univariate.csv' ))
    }

    if(grepl('feather', extension)){
      data_multivariate %>%
        write_feather(str_c(dir_output, '/', fn, '_event.feather' ))

      data_univariate %>%
        write_feather(str_c(dir_output, '/', fn, '_univariate.feather' ))
    }
  }
}

