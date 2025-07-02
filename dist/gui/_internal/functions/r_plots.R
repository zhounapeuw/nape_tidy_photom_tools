# Title     : f_general_functions.R
# Objective : collect functions related to plotting in R including general asethetics and specific plot generation
# Created by: Adam Gordon-Fennell (agg2248@uw.edu)
#             Garret Stuber Lab, University of Washington
# Created on: 09/02/2020

# themes ---------------------------------------------------------------------------------------------------------------
theme_ag01 <- function () {
  aes.axis_line_size <- 0.25
  aes.axis_font_size <- 8
  aes.axis_title_font_size <- 8
  aes.title_font_size <- 8

  theme_bw(base_size=10) %+replace%
    theme(
      axis.line = element_line(colour = "black", size = aes.axis_line_size), # set x / y axis
      axis.ticks = element_line(size = aes.axis_line_size),
      axis.ticks.length = unit(1, "mm"),
      axis.title = element_text(color = "black", size = aes.axis_title_font_size),
      axis.text.x = element_text(color = "black", size = aes.axis_font_size, hjust = 0.5),
      axis.text.y = element_text(color = "black", size = aes.axis_font_size, hjust = 1, vjust = 0.3, margin = margin(r = 1)),
      legend.title = element_text(color = "black", size = aes.axis_title_font_size),
      legend.text = element_text(color = "black", size = aes.axis_font_size),
      plot.title = element_text(color = "black", size = aes.title_font_size, hjust = 0.5),
      panel.grid.major = element_blank(), # turn off major grid
      panel.grid.minor = element_blank(), # turn off minor grid
      panel.border = element_blank(),     # turn of background panel border
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(color = 'black', size = aes.axis_title_font_size)
    )
}

theme_ag02 <- function () {
  aes.axis_line_size <- 0.25
  aes.axis_font_size <- 8
  aes.axis_title_font_size <- 8
  aes.title_font_size <- 8

  theme_bw(base_size=10) %+replace%
    theme(
      axis.line = element_line(colour = "black", size = aes.axis_line_size), # set x / y axis
      axis.ticks = element_line(size = aes.axis_line_size),
      axis.ticks.length = unit(1, "mm"),
      axis.title = element_text(color = "black", size = aes.axis_title_font_size),
      axis.text.x = element_text(color = "black", size = aes.axis_font_size, hjust = 0.5),
      axis.text.y = element_text(color = "black", size = aes.axis_font_size, hjust = 1, vjust = 0.3, margin = margin(r = 1)),
      legend.title = element_text(color = "black", size = aes.axis_title_font_size),
      legend.text = element_text(color = "black", size = aes.axis_font_size),
      plot.title = element_text(color = "black", size = aes.title_font_size, hjust = 0.5, face="bold"),
      panel.border = element_blank(),     # turn of background panel border
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(color = 'black', size = aes.axis_title_font_size)
    )
}


theme_ag_large <- function () {
  aes.axis_line_size <- 1
  aes.axis_font_size <- 18
  aes.axis_title_font_size <- 18
  aes.title_font_size <- 18

  theme_bw(base_size=18) %+replace%
    theme(
      axis.line = element_line(colour = "black", size = aes.axis_line_size), # set x / y axis
      axis.ticks = element_line(size = aes.axis_line_size),
      axis.ticks.length=unit(.25, "cm"),
      axis.title = element_text(color = "black", size = aes.axis_title_font_size),
      axis.text.x = element_text(color = "black", size = aes.axis_font_size, hjust = 0.5),
      axis.text.y = element_text(color = "black", size = aes.axis_font_size),
      legend.title = element_text(color = "black", size = aes.axis_title_font_size),
      legend.text = element_text(color = "black", size = aes.axis_font_size),
      plot.title = element_text(color = "black", size = aes.title_font_size, hjust = 0.5),
      panel.grid.major = element_blank(), # turn off major grid
      panel.grid.minor = element_blank(), # turn off minor grid
      panel.border = element_blank(),     # turn of background panel border
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(color = 'black', size = aes.axis_title_font_size)
    )
}

theme_ag_raster <- function () {
  aes.axis_line_size <- 0.25
  aes.axis_font_size <- 8
  aes.axis_title_font_size <- 8
  aes.title_font_size <- 8

  theme_bw(base_size=10) %+replace%
    theme(
      axis.ticks = element_line(size = aes.axis_line_size),
      axis.title = element_text(color = "black", size = aes.axis_title_font_size),
      legend.title = element_text(color = "black", size = aes.axis_title_font_size),
      legend.text = element_text(color = "black", size = aes.axis_font_size),
      plot.title = element_text(color = "black", size = aes.title_font_size, hjust = 0.5, face="bold"),
      panel.border = element_blank(),     # turn of background panel border
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(color = 'black', size = aes.axis_title_font_size),
      panel.grid.minor.y=element_blank(),
      panel.grid.major.y=element_blank(),
      axis.line.x = element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    )
}

# misc -----------------------------------------------------------------------------------------------------------------
remove_x_all <- function(plt){
  plt +
    theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

}

remove_y_all <- function(plt){
  plt +
    theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

}

return_n_x <- function(df, var_x){
  n_days <- df %>%
    select(blockname, subject, var_x) %>%
    unique() %>%
    group_by(subject) %>%
    summarise(day_count = n()) %>%
    pull(day_count) %>%
    max()

  return(n_days)
}

# geoms ----------------------------------------------------------------------------------------------------------------
# source: wilkelab/ungeviz
geom_hpline <- function(mapping = NULL, data = NULL,
                        stat = "identity", position = "identity",
                        ...,
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomHpline,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}

# source: wilkelab/ungeviz
GeomHpline <- ggproto("GeomHpline", GeomSegment,
  required_aes = c("x", "y"),
  non_missing_aes = c("size", "colour", "linetype", "width"),
  default_aes = aes(
    width = 0.5, colour = "black", size = 2, linetype = 1,
    alpha = NA
  ),

  draw_panel = function(self, data, panel_params, coord, arrow = NULL, arrow.fill = NULL,
                        lineend = "butt", linejoin = "round", na.rm = FALSE) {
    data <- mutate(data, x = x - width/2, xend = x + width, yend = y)
    ggproto_parent(GeomSegment, self)$draw_panel(
      data, panel_params, coord, arrow = arrow, arrow.fill = arrow.fill,
      lineend = lineend, linejoin = linejoin, na.rm = na.rm
    )
  }
)

# saving ---------------------------------------------------------------------------------------------------------------
save_pdf <- function(dir, fn, w, h){
  ggsave(
    filename =  str_c(dir, fn),
    device = NULL,
    path = NULL,
    scale = 1,
    width = w,
    height = h,
    units = c("in"),
    dpi = 300,
    useDingbats=FALSE
  )
}

# general functions ----------------------------------------------------------------------------------------------------
plt_string_facet <- function(plt, var_facet_y, var_facet_x, scales, facet_spacing){
  # returns plot faceted based on strings for y and x faceting variables
  # input NA in place of strings to avoid faceting

  if(!is.na(var_facet_y) & !is.na(var_facet_x)){
    plt <- plt + facet_grid(as.formula(str_c(var_facet_y, '~', var_facet_x)), scales = scales) +
      theme(panel.spacing = unit(facet_spacing, "lines"))
  }

  if(!is.na(var_facet_y) & is.na(var_facet_x)){
    plt <- plt + facet_grid(as.formula(str_c(var_facet_y, '~.')), scales = scales) +
      theme(panel.spacing = unit(facet_spacing, "lines"))
  }

  if(is.na(var_facet_y) & !is.na(var_facet_x)){
    plt <- plt + facet_grid(as.formula(str_c('.~', var_facet_x)), scales = scales) +
      theme(panel.spacing = unit(facet_spacing, "lines"))
  }

  return(plt)
}

plt_manual_dims <- function(plt, plt_dims){
  if(length(plt_dims) > 1){
    plt <- plt +
    force_panelsizes(rows = unit(plt_dims[1], "cm"),
                        cols = unit(plt_dims[2], "cm"))
  }

  return(plt)
}

plt_manual_scale_cartesian <- function(plt, plt_manual_scale_x, plt_manual_scale_y){
  # manually scale axis

  if(sum(is.na(plt_manual_scale_x) == 0) & sum(is.na(plt_manual_scale_y) == 0)){
    plt <- plt +
      coord_cartesian(xlim = plt_manual_scale_x, ylim = plt_manual_scale_y, expand = FALSE, clip = 'off')
  }

  if(sum(is.na(plt_manual_scale_x) == 0) & !sum(is.na(plt_manual_scale_y) == 0)){
    plt <- plt +
      coord_cartesian(xlim = plt_manual_scale_x, expand = FALSE, clip = 'off')
  }

  if(!sum(is.na(plt_manual_scale_x) == 0) & sum(is.na(plt_manual_scale_y) == 0)){
    plt <- plt +
      coord_cartesian(ylim = plt_manual_scale_y, expand = FALSE, clip = 'off')
  }

  return(plt)
}



# universal plots ----------------------------------------------------------------------------------------------------
plt_raster <- function(df, x, y, id, facet_x, facet_y, xlab, ylab, plt_y_lims, plt_x_lims, plt_dims, var_color, plt_scale_color){

  plt <- df %>%
    ungroup() %>%
    mutate(event_unique = row_number()) %>%
    select(na.omit(c(facet_y, facet_x, x, var_color, 'event_unique', y))) %>%
    mutate(y_dummy = !!as.name(y) + 1) %>%
    gather('id', 'y_pos', !!as.name(y):y_dummy) %>%

    ggplot(aes(event_ts_rel, y_pos, group = event_unique))

  if(is.na(var_color)) {plt <- plt + geom_line(size = 0.25)}
  if(!is.na(var_color)){plt <- plt + geom_line(size = 0.25, aes(color = !!as.name(var_color)))}


  if(!is.na(plt_scale_color) & !is.na(var_color)){
    if(plt_scale_color == 'mako'){plt <- plt + scale_color_viridis_d(option = 'mako', end = 0.9)}
    if(plt_scale_color == 'inferno'){plt <- plt + scale_color_viridis_d(option = 'inferno', end = 0.9)}
    if(plt_scale_color == 'viridis'){plt <- plt + scale_color_viridis_d(option = 'viridis', end = 0.9)}
  }

  if(!is.na(facet_x) & !is.na(facet_y)){plt <- plt + facet_grid(eval(expr(!!ensym(facet_y) ~ !!ensym(facet_x))))}
  if(!is.na(facet_x) &  is.na(facet_y)){plt <- plt + facet_grid(eval(expr(. ~ !!ensym(facet_x))))}
  if( is.na(facet_x) & !is.na(facet_y)){plt <- plt + facet_grid(eval(expr(!!ensym(facet_y) ~ .)))}

  plt <- plt +
    theme_ag01() +
    coord_cartesian(ylim = plt_y_lims, xlim = plt_x_lims, expand = FALSE) +
    force_panelsizes(rows = unit(plt_dims[1], "cm"),
                     cols = unit(plt_dims[2], "cm")) +
    xlab(xlab) +
    ylab(ylab)


  plt
}

plt_heatmap_trial_split <- function(df, var_x, var_fill, var_facet, var_trial, limits_fill, bin_seq, plt_dim, return_df){
  if(!hasArg(return_df)){
    return_df <- 0
  }

  # return heat plot of mean values over binned trials

  # setup split data
  df_split <- df %>%
    mutate(trial_split = cut(!!as.name(var_trial), bin_seq, label = bin_seq[1:length(bin_seq)-1])) %>%
    group_by_at(na.omit(c('trial_split', var_x, var_facet))) %>%
    summarise_at(vars(var_fill), mean) %>%
    mutate(trial_split = trial_split %>% as.character() %>% as.double()) %>%
    ungroup() %>%
    mutate(trial_split = trial_split + (bin_seq[2]/2)) %>%
    filter(!is.na(!!as.name(var_x)))

  # save solution count
  var_x_count <- df_split %>% select(!!as.name(var_x)) %>% unique() %>% nrow()

  plt <- df_split %>%
   ggplot(aes_string(x = var_x, y = 'trial_split', fill = var_fill)) +
   geom_tile() +
   theme_ag01() +
   coord_cartesian(ylim = c(0,max(bin_seq)), xlim = c(0.5, var_x_count + 0.5), expand = FALSE) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
   ylab('Trial Bin')

  # apply limits to plot
  if(sum(!is.na(limits_fill)) == 0){
   plt <- plt + scale_fill_continuous(low = 'black', high = 'white', oob = scales::squish)
  } else {
   plt <- plt + scale_fill_continuous(low = 'black', high = 'white', limits = limits_fill, oob = scales::squish)
  }

  # facet plot
  if(sum(!is.na(var_facet)) > 0){
   plt <- plt + facet_grid(as.formula(paste(".~", var_facet)))
  }
  plt <- plt %>%
    plt_manual_dims(plt_dim)

  #return plot
  if(return_df){
    return(list(plt, df_split))
  } else {
    return(plt)
  }
}


plt_mean_trial_split <- function(df, var_y, var_facet, var_trial, limits_fill, bin_seq, plt_dim, return_df){
  if(!hasArg(return_df)){
    return_df <- 0
  }

  # return heat plot of mean values over binned trials

  # setup split data
  df_split <- df %>%
    mutate(trial_split = cut(!!as.name(var_trial), bin_seq, label = bin_seq[1:length(bin_seq)-1])) %>%
    group_by_at(na.omit(c('trial_split', 'subject', var_facet))) %>%
    summarise_at(vars(var_y), mean) %>%
    mutate(trial_split = trial_split %>% as.character() %>% as.double()) %>%
    ungroup() %>%
    mutate(trial_split = trial_split + (bin_seq[2]/2))


  # save solution count

  plt <- df_split %>%
   ggplot(aes_string(x = 'trial_split', y = var_y)) +
   geom_line(alpha = 1/3, size = 0.25, aes(group = subject)) +
   stat_summary(fun = 'mean', geom = 'line', aes(group = 1), size = 0.25) +
   stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0, aes(group = 1), size = 0.25) +
   theme_ag01() +
   coord_cartesian(xlim = c(0,max(bin_seq)), expand = FALSE) +
   xlab('Trial Bin')

  # apply limits to plot
  if(sum(!is.na(limits_fill)) == 0){
   plt <- plt + scale_fill_continuous(low = 'black', high = 'white', oob = scales::squish)
  } else {
   plt <- plt + scale_fill_continuous(low = 'black', high = 'white', limits = limits_fill, oob = scales::squish)
  }

  # facet plot
  if(sum(!is.na(var_facet)) > 0){
   plt <- plt + facet_grid(as.formula(paste(".~", var_facet)))
  }
  plt <- plt %>%
    plt_manual_dims(plt_dim)

  #return plot
  if(return_df){
    return(list(plt, df_split))
  } else {
    return(plt)
  }

}

# fiber photometry -----------------------------------------------------------------------------------------------------
wavelength_to_rgb <- function(wavelength, gamma=0.8){

#    Original source:
#     https://gist.github.com/friendly/67a7df339aa999e2bcfcfec88311abfc#file-wavelength_to_rgb-r-L17
#    Based on code by Dan Bruton
#    http://www.physics.sfasu.edu/astro/color/spectra.html
#    '''
    if (wavelength >= 380 & wavelength <= 440) {
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ^ gamma
        G = 0.0
        B = (1.0 * attenuation) ^ gamma
        }
    else if (wavelength >= 440 & wavelength <= 490) {
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ^ gamma
        B = 1.0
        }
    else if (wavelength >= 490 & wavelength <= 510) {
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ^ gamma
        }
    else if (wavelength >= 510 & wavelength <= 580) {
        R = ((wavelength - 510) / (580 - 510)) ^ gamma
        G = 1.0
        B = 0.0
        }
    else if (wavelength >= 580 & wavelength <= 645) {
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ^ gamma
        B = 0.0
        }
    else if (wavelength >= 645 & wavelength <= 750) {
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ^ gamma
        G = 0.0
        B = 0.0
        }
    else {
        R = 0.0
        G = 0.0
        B = 0.0
        }
    R = R * 255
    G = G * 255
    B = B * 255
    return (rgb(floor(R), floor(G), floor(B), max=255))
}



plt_summary_number_of_files_per_subject <- function(df, plt_height, plt_width){
  # return point plot with the number of files for each subject

  plt <- df %>%
    select(blockname, subject, procedure) %>%
    unique() %>%
    group_by(subject, procedure) %>%
    summarise(count = n()) %>%
    ggplot(aes(subject, count, color = subject)) +
    geom_point() +
    scale_color_manual(values = polychrome(36) %>% as.vector()) +
    theme_ag01() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
    xlab('Subject') +
    ylab('Number of Files')


  if(!is.na(plt_height) & !is.na(plt_width)){
    plt <- plt +
      force_panelsizes(rows = unit(plt_height, "cm"),
                    cols = unit(plt_width, "cm"))
  }

  return(plt)
}

plt_summary_solution_in_data <- function(df, plt_height, plt_width){

plt <- df %>%
  select(blockname, subject, solution) %>%
  unique() %>%
  ggplot(aes(blockname, solution, fill = solution)) +
  geom_tile() +
  scale_fill_brewer(palette = 'Dark2') +
  theme_ag01() +
  theme(
    axis.text.x = element_text(angle = 90, hjust =0.5, vjust = 0.2),
    axis.title.x = element_blank(),
    axis.text.y = element_blank()
  ) +
  xlab('Blockname') +
  ylab('Solution') +
  ggtitle('Solutions in Data')

  if(!is.na(plt_height) & !is.na(plt_width)){
    plt <- plt +
      force_panelsizes(rows = unit(plt_height, "cm"),
                    cols = unit(plt_width, "cm"))
  }

return(plt)
}

plt_summary_signal_ids_in_data <- function(df, plt_height, plt_width){
  plt_data <- df %>%
    select(blockname, fiber_id, signal_wavelength, control_wavelength) %>%
    unique() %>%
    group_by(blockname) %>%
    mutate(signal_number = row_number()) %>%
    gather('id', 'wavelength', signal_wavelength:control_wavelength) %>%
    arrange(blockname, wavelength)

  plt_colors <- plt_data %>%
    select(wavelength) %>%
    unique() %>%
    rowwise() %>%
    mutate(wavelength = wavelength_to_rgb(wavelength)) %>%
    pull(wavelength)

  plt <- plt_data %>%
    mutate(wavelength = as.character(wavelength),
           fiber_id = str_c('fiber_', fiber_id),
           signal_number = str_c('signal_', signal_number),
           id = id %>% str_remove('_wavelength')
           ) %>%
    ggplot(aes(blockname, id, fill = wavelength)) +
      geom_tile() +
      scale_fill_manual(values = plt_colors) +
      #scale_fill_brewer(palette = 'Accent') +
      facet_grid(interaction(fiber_id, signal_number)~.) +
      theme_ag01() +
      theme(axis.text.x = element_text(angle = 90, hjust =0.5, vjust = 0.2),
            axis.title.x = element_blank()) +
      ylab('Siginal IDs') +
      ggtitle('Signals in Data')

  if(!is.na(plt_height) & !is.na(plt_width)){
    plt <- plt +
      force_panelsizes(rows = unit(plt_height, "cm"),
                    cols = unit(plt_width, "cm"))
  }

return(plt)
}

plt_summary_licking_behavior <- function(df, plt_height, plt_width){
  plt_data <- df %>%
    select(blockname, subject, solution, trial_num, lick_count) %>%
    unique() %>%
    group_by(subject, solution) %>%
    summarise(lick_count_mean = lick_count %>% mean(),
              lick_count_proportion = sum(lick_count > 0) / n())


  plt01 <- plt_data %>%
    ggplot(aes(solution, lick_count_mean, color = solution, group = 1)) +
    geom_point(size = 3, shape = 21, stroke = 0.25, fill = NA, alpha = 0.5) +
    stat_summary(fun = 'mean', geom = 'hpline', size = 0.5) +
    stat_summary(fun.data = 'mean_se', geom = 'errorbar', size = 0.5, width = 0) +
    theme_ag01() +
    scale_color_brewer(palette = 'Dark2') +
    theme(
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    ylab('Mean Lick Count')


  plt02 <- plt_data %>%
    ggplot(aes(solution, lick_count_proportion, color = solution, group = 1)) +
    geom_point(size = 3, shape = 21, stroke = 0.25, fill = NA, alpha = 0.5) +
    stat_summary(fun = 'mean', geom = 'hpline', size = 0.5) +
    stat_summary(fun.data = 'mean_se', geom = 'errorbar', size = 0.5, width = 0) +
    theme_ag01() +
    scale_color_brewer(palette = 'Dark2') +
    theme(
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    ylab('Trial Lick Proportion')

  plt03 <- plt_data %>%
    ggplot(aes(solution, subject, fill = lick_count_mean)) +
    geom_tile() +
    scale_fill_continuous(low = 'black', high = 'white', oob = scales::squish) +
    theme_ag01() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.2)) +
    xlab('Solution') +
    ylab('Subject')

  if(!is.na(plt_height) & !is.na(plt_width)){
    plt01 <- plt01 +
      force_panelsizes(rows = unit(plt_height, "cm"),
                    cols = unit(plt_width, "cm"))

    plt02 <- plt02 +
      force_panelsizes(rows = unit(plt_height, "cm"),
                    cols = unit(plt_width, "cm"))

    plt03 <- plt03 +
      force_panelsizes(rows = unit(plt_height, "cm"),
                    cols = unit(plt_width, "cm"))

  }

return(list(plt01, plt02, plt03))
}

plt_combined_summary_multi_spout <- function(df, plt_height, plt_width){
  df <- df %>%
    group_by(blockname, fiber_id, signal_wavelength, trial_num) %>%
    filter(time_rel == min(time_rel)) %>%
    ungroup()

  plt01 <- df %>%
    plt_summary_signal_ids_in_data(NA,NA)

  plt02 <- df %>%
    plt_summary_solution_in_data(NA,NA)

  plt03 <- df %>%
    plt_summary_number_of_files_per_subject(NA,NA)

  plt04_list <- df %>%
    plt_summary_licking_behavior(NA, NA)

  plt_04_01 <- plt04_list[[1]]

  plt_04_02 <- plt04_list[[2]]

  plt_04_03 <- plt04_list[[3]]

  num_plt <- 6

  plt01 + plt02 + plt03 + plt_04_01 + plt_04_02 + plt_04_03 +
    plot_layout(
      ncol = 1,
      heights = unit(rep(plt_height, num_plt), rep('cm', num_plt)),
      widths =  unit(c(plt_width), c('cm')),
      )


}


plt_streams_peth_multispout_trace <- function(df, y, viridis_scale, var_facet_y, var_facet_x, plt_dims){
  plt <- df %>%
    ggplot(aes_string('time_rel', y, color = 'solution', fill = 'solution', group = 'solution')) +
    stat_summary(fun.data = 'mean_se', geom = 'ribbon', color = NA, alpha = 1/3) +
    stat_summary(fun = 'mean', geom = 'line') +
    geom_vline(xintercept = c(0,3), linetype = 2, alpha = 0.3) +
    scale_color_viridis_d(option = viridis_scale, end = 0.9) +
    scale_fill_viridis_d(option = viridis_scale, end = 0.9) +
    theme_ag01() +
    xlab('Time (s)')


  plt <- plt %>%
    plt_string_facet(var_facet_y, var_facet_x, 'free', 1)

  plt <- plt %>%
    plt_manual_dims(plt_dims)



  plt
}



plt_fp_multispout_heatplts_combined <- function(streams_peth, data_trial_summary, fp_experiment, fp_procedure, log_fp){
  import_fns <- log_fp %>%
    filter(experiment == fp_experiment) %>%
    filter(procedure == fp_procedure) %>%
    arrange(blockname) %>%
    pull(blockname)

  plt_licks <- plt_heatmap_trial_split(
    df = data_trial_summary %>% filter(file_name %in% import_fns),
    var_x = 'solution',
    var_fill = 'lick_count',
    var_facet = NA,
    var_trial = 'trial_num',
    limits_fill = c(0,16),
    bin_seq = seq(0, 100, 25),
    plt_dim = NA)

  plt_licks <- plt_licks +
    ggtitle('Trial Licks')

  df <- streams_peth %>%
    filter(time_rel > 0, time_rel < 3) %>%
    group_by(subject, group, channel_id, solution, trial_num) %>%
    summarise(delta_signal_poly_zscore = delta_signal_poly_zscore %>% mean()) %>%

    filter(!(group == 'green' & str_detect(channel_id, '560'))) %>%
    filter(!(group == 'red'   & str_detect(channel_id, '465')))


  plt_fluor <- plt_heatmap_trial_split(
    df = df,
    var_x = 'solution',
    var_fill = 'delta_signal_poly_zscore',
    var_facet = 'channel_id',
    var_trial = 'trial_num',
    limits_fill = NA,
    bin_seq = seq(0, 100, 25),
    plt_dim = NA)

  plt_fluor <- plt_fluor +
    ggtitle('Trial Mean Fluor')

  return(list(plt_licks, plt_fluor))
}

plt_fp_multispout_mean_trialsplit_combined <- function(streams_peth, data_trial_summary, fp_experiment, fp_procedure, log_fp){
  import_fns <- log_fp %>%
    filter(experiment == fp_experiment) %>%
    filter(procedure == fp_procedure) %>%
    arrange(blockname) %>%
    pull(blockname)

  plt_licks <- plt_mean_trial_split(
    df = data_trial_summary %>% filter(file_name %in% import_fns),
    var_y = 'lick_count',
    var_facet = NA,
    var_trial = 'trial_num',
    limits_fill = c(0,16),
    bin_seq = seq(0, 100, 25),
    plt_dim = NA)

  plt_licks <- plt_licks +
    ggtitle('Trial Licks')

  df <- streams_peth %>%
    filter(time_rel > 0, time_rel < 3) %>%
    group_by(subject, group, channel_id, solution, trial_num) %>%
    summarise(delta_signal_poly_zscore = delta_signal_poly_zscore %>% mean()) %>%

    filter(!(group == 'green' & str_detect(channel_id, '560'))) %>%
    filter(!(group == 'red'   & str_detect(channel_id, '465')))


  plt_fluor <- plt_mean_trial_split(
    df = df,
    var_y = 'delta_signal_poly_zscore',
    var_facet = 'channel_id',
    var_trial = 'trial_num',
    limits_fill = NA,
    bin_seq = seq(0, 100, 25),
    plt_dim = NA)

  plt_fluor <- plt_fluor +
    ggtitle('Trial Mean Fluor')


  return(list(plt_licks, plt_fluor))

}


# head-fixed operant ---------------------------------------------------------------------------------------------------
plt_rotary_training_criteria <- function(df_rotary_summary, var_x, var_y_set, var_group, color_values, n_day_summary){

  if(var_y_set == 'criteria_counts'){
    data_plt <- df_rotary_summary %>%
      filter(event_id_char %in% c('active_rotation_criteria', 'inactive_rotation_criteria')) %>%
      mutate(event_id_char = ifelse(event_id_char == 'active_rotation_criteria', 'active', 'inactive'))

    # extract max y value for aesthetics
    y_max <- data_plt %>%
      pull(event_count) %>%
      max()

    var_y <- 'event_count'
  }

  if(var_y_set == 'rotation_turn'){
    data_plt <- df_rotary_summary %>%
      filter(event_id_char %in% c('active_rotation', 'inactive_rotation')) %>%
      mutate(event_id_char = ifelse(event_id_char == 'active_rotation', 'active', 'inactive'))

    # extract max y value for aesthetics
    y_max <- data_plt %>%
      pull(rotation_turn) %>%
      max()

    var_y <- 'rotation_turn'
  }



  # number of days for aesthetics
  n_x <- return_n_x(data_plt, var_x)


  if(!is.na(var_group)){
    plt_rotary_training_criteria <- data_plt %>%
      ggplot(aes_string(var_x, var_y, group = interaction('event_id_char', var_group), linetype = 'event_id_char', color = var_group)) +
      geom_line(aes(group = interaction(subject, event_id_char)), alpha = 1/3) +
      stat_summary(fun = 'mean', geom = 'line',
                   aes(group = interaction(!!as.name(var_group), !!as.name('event_id_char')))) +
      stat_summary(fun.data = 'mean_se', geom = 'errorbar',
                   aes(group = interaction(!!as.name(var_group), !!as.name('event_id_char'))),
                   width = 0, linetype = 1) +
      scale_color_manual(values = color_values)
  } else {
    plt_rotary_training_criteria <- data_plt %>%
      ggplot(aes_string(var_x, var_y, group = 'event_id_char', linetype = 'event_id_char')) +
      geom_line(aes(group = interaction(subject, event_id_char)), alpha = 1/3) +
      stat_summary(fun = 'mean', geom = 'line',
                   aes(group = event_id_char)) +
      stat_summary(fun.data = 'mean_se', geom = 'errorbar',
                   aes(group = event_id_char),
                   width = 0, linetype = 1)
  }

  plt_rotary_training_criteria <- plt_rotary_training_criteria +
    theme_ag01() +
    scale_linetype_manual(values = c(1,2)) +
    coord_cartesian(ylim = c(0,y_max*1.1), xlim = c(0.5, n_x + 0.5), expand = FALSE) +
    xlab('Day')

  if(var_y_set == 'criteria_counts'){
    plt_rotary_training_criteria <- plt_rotary_training_criteria +
      ylab('Criteria Count') +
      ggtitle('Criteria Counts \n')
  }

  if(var_y_set == 'rotation_turn'){
    plt_rotary_training_criteria <- plt_rotary_training_criteria +
      ylab('Rotation (Turn)') +
      ggtitle('Rotation \n')
  }

  # plot n day mean
    data_plt_summary <- data_plt %>%
      mutate(day = day %>% str_remove('d') %>% as.integer()) %>%
      group_by(subject) %>%
      mutate(day = day - max(day)) %>%
      filter(day > -n_day_summary)

      if(is.na(var_group)){
        data_plt_summary <- data_plt_summary %>%
          group_by(subject, event_id_char) %>%
          summarise(var_y_mean = !!as.name(var_y) %>% mean())
      } else {
        data_plt_summary <- data_plt_summary %>%
          group_by(subject, event_id_char, !!as.name(var_group)) %>%
          summarise(var_y_mean = !!as.name(var_y) %>% mean()) %>%
          ungroup()

        n_groups <- data_plt_summary %>%
          select(!!as.name(var_group)) %>%
          unique() %>%
          nrow()
      }

    if(is.na(var_group)){
      plt_rotary_training_criteria_summary <- data_plt_summary %>%
        ggplot(aes(event_id_char, var_y_mean, group = event_id_char, linetype = event_id_char)) +
        geom_quasirandom(width = 0.1, shape = 21, fill = NA, size = 3, stroke = 0.25, alpha = 1/3) +
        stat_summary(fun = 'mean', geom = 'hpline', size = 0.5) +
        stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0, size = 0.5, linetype = 1) +
        coord_cartesian(ylim = c(0,y_max*1.1), xlim = c(0, 3), expand = FALSE, clip = 'off')
        } else {
      plt_rotary_training_criteria_summary <- data_plt_summary %>%
        ggplot(aes(event_id_char, var_y_mean, group = !!as.name(var_group), color = !!as.name(var_group), linetype = 'event_id_char')) +
        geom_quasirandom(width = 0.1, shape = 21, fill = NA, size = 3, stroke = 0.25, alpha = 1/3) +
        stat_summary(fun = 'mean', geom = 'hpline', size = 0.5) +
        stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0, size = 0.5, linetype = 1) +
        scale_color_manual(values = color_values) +
        coord_cartesian(ylim = c(0,y_max*1.1), xlim = c(0, 3), expand = FALSE, clip = 'off') +
        facet_grid(as.formula(str_c('.~',var_group)))
      }

      plt_rotary_training_criteria_summary <- plt_rotary_training_criteria_summary +
        theme(legend.position = "none") +
        ggtitle(str_c(n_day_summary, 'd Mean\n')) +
        theme_ag01() %>%
        remove_x_all() %>%
        remove_y_all()




  return(list(plt_rotary_training_criteria, plt_rotary_training_criteria_summary))

}

plt_rotary_training_criteria_bias <- function(df_rotary_summary, var_x, var_y_set, var_group, color_values, n_day_summary){

  if(var_y_set == 'criteria_counts'){
    data_plt <- df_rotary_summary %>%
      filter(event_id_char %in% c('active_rotation_criteria', 'inactive_rotation_criteria')) %>%
      select(na.omit(c('blockname', 'subject', var_x, var_group, 'event_id_char', 'event_count'))) %>%
      spread(event_id_char, event_count) %>%
      mutate(inactive_rotation_criteria = ifelse(is.na(inactive_rotation_criteria), 0, inactive_rotation_criteria)) %>%
      mutate(active_bias = active_rotation_criteria / (active_rotation_criteria + inactive_rotation_criteria))
  }

  if(var_y_set == 'rotation_turn'){
    data_plt <- df_rotary_summary %>%
      filter(event_id_char %in% c('active_rotation', 'inactive_rotation')) %>%
      select(na.omit(c('blockname', 'subject', var_x, var_group, 'event_id_char', 'rotation_turn'))) %>%
      spread(event_id_char, rotation_turn) %>%
      mutate(active_bias = active_rotation / (active_rotation + inactive_rotation))
  }

  n_x <- return_n_x(data_plt, var_x)


  if(is.na(var_group)){
  plt_rotary_training_criteria_bias <- data_plt %>%
    ggplot(aes_string(var_x, 'active_bias', group = 1)) +
    geom_hline(yintercept = 0.5, color = 'grey') +
    geom_line(aes(group = subject), size = 0.25, alpha = 1/3) +
    stat_summary(fun = 'mean', geom = 'line', size = 0.5) +
    stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0, size = 0.5)
  } else {
  plt_rotary_training_criteria_bias <- data_plt %>%
    ggplot(aes_string(var_x, 'active_bias', group = var_group, color = var_group)) +
    geom_hline(yintercept = 0.5, color = 'grey') +
    geom_line(aes(group = subject), size = 0.25, alpha = 1/3) +
    stat_summary(fun = 'mean', geom = 'line', size = 0.5) +
    stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0, size = 0.5) +
    scale_color_manual(values = color_values)
  }

  plt_rotary_training_criteria_bias <- plt_rotary_training_criteria_bias +
    theme_ag01() +
    coord_cartesian(ylim = c(0,1), xlim = c(0.5, n_x + 0.5), expand = FALSE) +
    xlab('Day')


    if(var_y_set == 'criteria_counts'){
      plt_rotary_training_criteria_bias <- plt_rotary_training_criteria_bias +
        ylab('Active Criteria Bias') +
        ggtitle('Active Criteria Bias\n')
    }

    if(var_y_set == 'rotation_turn'){
      plt_rotary_training_criteria_bias <- plt_rotary_training_criteria_bias +
        ylab('Active Rotation Bias') +
        ggtitle('Active Rotation Bias\n')
    }

  # plot n day mean
    data_plt_summary <- data_plt %>%
      mutate(day = day %>% str_remove('d') %>% as.integer()) %>%
      group_by(subject) %>%
      mutate(day = day - max(day)) %>%
      filter(day > -n_day_summary)

      if(is.na(var_group)){
        data_plt_summary <- data_plt_summary %>%
          group_by(subject) %>%
          summarise(active_bias_mean = active_bias %>% mean())
      } else {
        data_plt_summary <- data_plt_summary %>%
          group_by(subject, !!as.name(var_group)) %>%
          summarise(active_bias_mean = active_bias %>% mean()) %>%
          ungroup()

        n_groups <- data_plt_summary %>%
          select(!!as.name(var_group)) %>%
          unique() %>%
          nrow()
      }

    if(is.na(var_group)){
      plt_rotary_training_criteria_bias_summary <- data_plt_summary %>%
        ggplot(aes(1, active_bias_mean, group = 1)) +
        geom_hline(yintercept = 0.5, color = 'grey') +
        geom_quasirandom(width = 0.1, shape = 21, fill = NA, size = 3, stroke = 0.25, alpha = 1/3) +
        stat_summary(fun = 'mean', geom = 'hpline', size = 0.5) +
        stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0, size = 0.5) +
        coord_cartesian(ylim = c(0,1), xlim = c(0, 2), expand = FALSE, clip = 'off')
      } else {
      plt_rotary_training_criteria_bias_summary <- data_plt_summary %>%
        ggplot(aes_string(var_group, 'active_bias_mean', group = var_group, color = var_group)) +
        geom_hline(yintercept = 0.5, color = 'grey') +
        geom_quasirandom(width = 0.1, shape = 21, fill = NA, size = 3, stroke = 0.25, alpha = 1/3) +
        stat_summary(fun = 'mean', geom = 'hpline', size = 0.5) +
        stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0, size = 0.5) +
        scale_color_manual(values = color_values) +
        coord_cartesian(ylim = c(0,1), xlim = c(0, n_groups + 1), expand = FALSE, clip = 'off')
      }

      plt_rotary_training_criteria_bias_summary <- plt_rotary_training_criteria_bias_summary +
        theme(legend.position = "none") +
        ggtitle(str_c(n_day_summary, 'd Mean\n')) +
        theme_ag01() %>%
        remove_x_all() %>%
        remove_y_all()


        if(var_y_set == 'criteria_counts'){
          plt_rotary_training_criteria_bias_summary <- plt_rotary_training_criteria_bias_summary +
            ylab('Active Criteria Bias')
        }

        if(var_y_set == 'rotation_turn'){
          plt_rotary_training_criteria_bias_summary <- plt_rotary_training_criteria_bias_summary +
            ylab('Active Rotation Bias') +
            ggtitle('Active Rotation Bias\n')
        }

  return(list(plt_rotary_training_criteria_bias, plt_rotary_training_criteria_bias_summary))

}

# head-fixed multi-spout brief access ----------------------------------------------------------------------------------

plt_lick_mean_trial_lick_count_ring <- function(df, var_x, var_y, var_facet_y, var_facet_x, facet_scales, facet_spacing, scale_viridis_option, plt_manual_scale_y, plt_manual_scale_x, plt_manual_dims){

  plt <- df %>%
    ggplot(aes_string(var_x, var_y, color = var_x)) +
    geom_quasirandom(size = 2, shape = 21, fill = NA, alpha = 0.5, width = 0.3) +
    stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0) +
    stat_summary(fun = 'mean', geom = 'hpline', aes(group = 1), size = 0.5, width = 0.5) +

    ylab('Mean Trial Lick Count') +

    scale_color_viridis_d(option = scale_viridis_option, end = 0.9) +
    theme_ag01() +
    theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())


  plt <- plt %>%
    plt_string_facet(var_facet_y, var_facet_x, facet_scales, facet_spacing) %>%
    plt_manual_scale_cartesian(plt_manual_scale_x, plt_manual_scale_y) %>%
    plt_manual_dims(plt_manual_dims)


  return(plt)

}

plt_lick_mean_trial_lick_count_line <- function(df, var_x, var_y, var_facet_y, var_facet_x, facet_scales, facet_spacing, plt_manual_scale_y, plt_manual_scale_x, plt_manual_dims){

  plt <- df %>%
    ggplot(aes_string(var_x, var_y)) +
    geom_line(aes(group = subject), size = 1/4, alpha = 1/4) +
    stat_summary(fun.data = 'mean_se', geom = 'errorbar', aes(group = 1), width = 0) +
    stat_summary(fun = 'mean', geom = 'line', aes(group = 1), size = 0.25) +

    ylab('Mean Trial Lick Count') +
    xlab('Solution') +

    theme_ag01() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1))


  plt <- plt %>%
    plt_string_facet(var_facet_y, var_facet_x, facet_scales, facet_spacing) %>%
    plt_manual_scale_cartesian(plt_manual_scale_x, plt_manual_scale_y) %>%
    plt_manual_dims(plt_manual_dims)


  return(plt)

}

plt_multi_spout_lick_summary_training <- function(data_spout_summary, limits_licks){

plt1 <- data_spout_summary %>%
  ggplot(aes(solution, lick_count_trial_mean)) +
  geom_line(aes(group = subject), alpha = 1/3)+
  stat_summary(fun = 'mean', geom = 'line', aes(group = 1)) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', aes(group = 1), width = 0) +
  facet_grid(.~day) +
  theme_ag01() +
  ylab('Mean Trial Lick Count')

plt1 <- plt1 %>%
  remove_x_all()


plt1_mean <-
  data_spout_summary %>%
  group_by(subject, solution) %>%
  summarise(lick_count_trial_mean = lick_count_trial_mean %>% mean()) %>%
  plt_lick_mean_trial_lick_count_line(
    df = .,
    var_x = 'solution',
    var_y = 'lick_count_trial_mean',
    var_facet_y = NA,
    var_facet_x = NA,
    facet_scales = 'free',
    facet_spacing = 1,

    plt_manual_scale_y = c(0,15),
    plt_manual_scale_x = c(0,6),
    plt_manual_dims = c(3,2)
  )

plt1_mean <- plt1_mean %>%
  remove_x_all() %>%
  remove_y_all() +
  ggtitle('Mean')


plt2 <- data_spout_summary %>%
  ggplot(aes(solution, subject, fill = lick_count_trial_mean)) +
  geom_tile()+
  facet_grid(.~day) +
  theme_ag01() +
  scale_fill_continuous(low = 'black', high = 'white', limits = limits_licks, oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  ylab('Subject')

plt2_mean <- data_spout_summary %>%
  group_by(subject, solution) %>%
  summarise(lick_count_trial_mean = lick_count_trial_mean %>% mean()) %>%
  ggplot(aes(solution, subject, fill = lick_count_trial_mean)) +
  geom_tile()+
  theme_ag01() +
  scale_fill_continuous(low = 'black', high = 'white', limits = limits_licks, oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.3)) +
  ylab('Subject')

plt2_mean <- plt2_mean %>%
  remove_y_all()

  n_day <- data_spout_summary %>%
    select(day) %>% unique() %>%
    nrow()

plt_combined <- plt1 + plt1_mean + plt2 + plt2_mean +
  plot_layout(ncol = 2,
              widths  = unit(c(n_day * 2.5,2.5), rep('cm', 2)),
              heights = unit(rep(3, 2), rep('cm', 2)),
              guides = 'collect'
              )
return(plt_combined)
}











