# Title     : r_rotary.R
# Objective : to house functions used for extraction / pre processing rotary files
# Created by: adam gordon-fennell
# Created on: 04/16/2020

# required functions
require(tidyverse)
require(lubridate)
require(readxl)
require(janitor)
require(feather)

extract_etho <- function(input_dirs, folder_output, save_trial_files, save_combined_file, file_name_subject_day, user_parameters, verbose){
  # extracts data and parameters saved within raw ethovision data file
  # output file types
  #   - trial parameters (1 row of parameters for the trial)
  #   - trial raw (time serises data)
  #   - combined parameters
  #   - combined raw
  # all files exported as .feather files in tidy format
  #
  # function parameters
  #   - input_dirs: vector of input directories to be extracted
  #   - folder_output: directory of outputs
  #   - save_trial_files: toggle for saving individual trials
  #   -  save_combined_file: toggle for saving combined dataframe
  #   - file_name_subject_day: toggle to name file using 'subject' and 'day' from pramaters; FALSE: name file using date and trial
  #   - user_parameters: vector of parameters within ethovision data to keep track of
  #   - verbose: toggle to print each file extracted


  n_combined_initial <- 1

  print('extracting files from following directories:')
  for(input_dir in input_dirs){
    print(input_dir)
  }

  print(' ')

  for(input_dir in input_dirs){

    file_names <- list.files(input_dir) %>%
      str_c(input_dir, '/', .)

    for(file_name in file_names){

    header_rows <- read_csv(file_name, n_max = 1, col_names = FALSE, col_types = cols()) %>% pull(X2)

    # read parameters
    parameters <- read_csv(file_name, n_max = header_rows-3, col_names = c('var', 'value', 'empty'), col_types = cols()) %>%
      select(-empty) %>%
      spread(var, value) %>%
      clean_names() %>%
        separate(video_start_time, into = c('date', 'video_start_time'), sep = ' ')

    if(parameters %>% pull(subject) != 'NA'){
      # read in raw data
      raw <- read_csv(file_name, skip = header_rows-1, na = "-", col_types = cols())

      raw <- raw %>%
        select(-colnames(raw)[dim(raw)[2]])

      raw_colnames <- read_csv(file_name, skip = header_rows-2, n_max = 1, col_types = cols()) %>%
        clean_names()

      raw_colnames <- raw_colnames %>%
        select(-colnames(raw_colnames)[dim(raw_colnames)[2]])

      colnames(raw) <- colnames(raw_colnames)

      # save categorical data from parameters

      if(sum(is.na(user_parameters))== 0){
        raw <- bind_cols(
          raw,
          parameters %>% select(date, trial_name, user_parameters) %>%
          mutate(count = c(dim(raw)[1])) %>%
          uncount(count)
          )
      }else{
          raw <- bind_cols(
          raw,
          parameters %>% select(date, trial_name) %>%
          mutate(count = c(dim(raw)[1])) %>%
          uncount(count)
          )
      }


      if(save_trial_files & dim(raw)[1]>1){
        if(file_name_subject_day){
          output_fn_trial <-  str_c(parameters %>% pull(day),
                                    '_',
                                    parameters %>% pull(subject))

        } else {
          output_fn_trial <- str_c(parameters %>% pull(date) %>% str_replace_all('/','_'),
                                   '_',
                                   parameters %>% pull(trial_name) %>% str_remove_all(' '))
        }
        if(verbose){print(str_c('saving trial file: ', output_fn_trial))}
        raw        %>% write_feather(str_c(folder_output, output_fn_trial, '_raw.feather'))
        parameters %>% write_feather(str_c(folder_output, output_fn_trial, '_parameters.feather'))
      }

      if(save_combined_file & dim(raw)[1]>1){
        if(n_combined_initial){
          raw_combined <- raw
          parameters_combined <- parameters
          n_combined_initial <- 0
        } else {
          raw_combined <- raw_combined %>% bind_rows(raw)
          parameters_combined <- parameters_combined %>% bind_rows(parameters)
        }

      if(dim(raw)[1]==1){
        print(str_c(
          'ERROR: MISSING TRACK INFO FOR: ',
          parameters %>% pull(date) %>% str_replace_all('/','_'),
          ' ',
          parameters %>% pull(trial_name) %>% str_remove_all(' '))
        )
      }

      }
    } # subject != NA
    } # file loop
  } # dir loop

  if(save_combined_file){
   if(verbose){print('saving combined raw and combined parameters')}
    raw_combined        %>% write_feather(str_c(folder_output, 'combined_raw.feather'))
    parameters_combined %>% write_feather(str_c(folder_output, 'combined_parameters.feather'))

  }

    print('extraction complete')
}
