# Title     : f_general_functions.R
# Objective : collect functions for performing stats and tidying oputputs
# Created by: Adam Gordon-Fennell (agg2248@uw.edu)
#             Garret Stuber Lab, University of Washington
# Created on: 09/10/2020

require(broom)

####### general functions
# generate labels for significance
quick_sig <- function(stats_df){
    stats_df %>%
      mutate(sig_text = ifelse(p_value < 0.05, 'sig', ''),
             sig_symbol = ifelse(p_value < 0.001, '***',
                                 ifelse(p_value < 0.01, '**',
                                        ifelse(p_value < 0.05, '*',
                                               'n.s.'))))
}

# save aov
save_aov <- function(dir_output, prefix, anova_summary, pairs_hsd, df){
  if(!is.na(dir_output)){
    fn <- str_c(dir_output, '/', prefix, '_aov.csv')
    anova_summary %>% write_csv(fn)
    print(str_c('saved file: ', fn))

    fn <- str_c(dir_output, '/', prefix, '_aov_hsd.csv')
    pairs_hsd %>% write_csv(fn)
    print(str_c('saved file: ', fn))
    
    fn <- str_c(dir_output, '/', prefix, '_data.csv')
    df %>% write_csv(fn)
    print(str_c('saved file: ', fn))
  }
}


save_ttest <- function(dir_output, prefix, stats_result){
  if(!is.na(dir_output)){
  fn <- str_c(dir_output, '/', prefix, '_ttest.csv')
  stats_result %>% write_csv(fn)
  print(str_c('saved file: ', fn))
  }
 }

####### specific stats functions
# ks test
do_ks <- function(df, var1, var2){
    ks.test(
            df %>% pull(var1),
            df %>% pull(var2),
            alternative = 'two.sided'
    ) %>%
    broom::tidy() %>%
    mutate(ag_stat = 'do_ks')
}

tidy_ks <- function(df, var1, var2, ...){
    df %>%
      group_by_(...) %>%
      do(do_ks(., var1, var2))
}

# ttest paired
do_ttest_paired <- function(df, var1, var2){
    t.test(
            df %>% pull(var1),
            df %>% pull(var2),
            paired = TRUE,
            alternative = 'two.sided'
    ) %>%
    broom::tidy() %>%
    mutate(ag_stat = 'do_ttest_paired')
}

tidy_ttest_paired <- function(df, dir_output, prefix, var1, var2, ...){
    stats_result <- df %>%
      group_by_(...) %>%
      do(do_ttest_paired(., var1, var2)) %>%
      mutate(var1 = var1,
             var2 = var2) %>%
      select(var1, var2, everything()) %>%
      clean_names()

  save_ttest(dir_output, prefix, stats_result)

  return(stats_result)
}

# ttest unpaired
do_ttest_unpaired <- function(df, var1, var2){
    t.test(
            df %>% pull(var1),
            df %>% pull(var2),
            paired = FALSE,
            alternative = 'two.sided'
    ) %>%
    broom::tidy() %>%
    mutate(ag_stat = 'do_ttest_unpaired')
}

tidy_ttest_unpaired <- function(df, dir_output, prefix, var1, var2, ...){
    stats_result <- df %>%
      group_by_(...) %>%
      do(do_ttest_unpaired(., var1, var2)) %>%
      mutate(var1 = var1,
             var2 = var2) %>%
      select(..., var1, var2, everything()) %>%
      clean_names()

  save_ttest(dir_output, prefix, stats_result)

  return(stats_result)
}

do_ttest_unpaired_grouped <- function(df, var_iv, var_dv){
  # similar to do_ttest_unpaired, but takes in the column name for the independent / dependent variables instead

    df_var_iv_levels <- df %>%
      pull(var_iv) %>%
      unique()

    df <- df %>%
      spread(!!as.name(var_iv), !!as.name(var_dv))

    t.test(
            df %>% pull(!!as.name(df_var_iv_levels[1])),
            df %>% pull(!!as.name(df_var_iv_levels[2])),
            paired = FALSE,
            alternative = 'two.sided'
    ) %>%
    broom::tidy() %>%
    clean_names() %>%
    mutate(ag_stat = 'do_ttest_unpaired',
           var1    = df_var_iv_levels[1],
           var2    = df_var_iv_levels[2]) %>%
    select(ag_stat, var1, var2, everything())
}


tidy_ttest_unpaired_grouped <- function(df, dir_output, prefix, var_dependent, value_dependent, var_between, value_between, ...){
  # performs ttest on grouped data without having to include variable names for the independent variables you wish to compare

  check_count_value_between <- df %>%
    ungroup() %>%
    select(var_between, var_dependent, value_between, ...) %>%
    unique() %>%
    group_by(var_between, var_dependent, ...) %>%
    summarise(count_value_between = n()) %>%
    filter(count_value_between != 2)

  if(nrow(check_count_value_between) > 0){
    print('tidy_ttest_unpaired_grouped stopped... Some independent variables include more than 2 levels')
    print(check_count_value_between)
    break
  }

  stats_result <- df %>%
    group_by_(var_dependent, var_between, ...) %>%
    do(do_ttest_unpaired_grouped(., value_between, value_dependent)) %>%
    quick_sig()


  save_ttest(dir_output, prefix, stats_result)

  return(stats_result)
}


# pearson correlation
do_cor <- function(df, var1, var2){
    tidy_cor_test <- cor.test(
            df %>% pull(var1),
            df %>% pull(var2),
            method = 'pearson'
    ) %>%
    broom::tidy() %>%
    mutate(ag_stat = 'do_cor') %>%
      rename('p_value' = p.value,
             'conf_low' = conf.low,
             'con_high' = conf.high)

    model <- lm(df %>% pull(var2) ~ df %>% pull(var1))

    model %>%
      broom::tidy() %>%
      clean_names() %>%
      mutate(parameter = ifelse(term == '(Intercept)', 'intercept', 'slope')) %>%
      mutate(parameter = str_c('estimate_', parameter)) %>%
      select(parameter, estimate) %>%
      spread(parameter, estimate) %>%
      bind_cols(tidy_cor_test)

}

tidy_cor <- function(df, var1, var2, ...){
    df %>%
      group_by_(...) %>%
      do(do_cor(., var1, var2)) %>%
      quick_sig()
}

# cross correlation
do_cross_correlation <- function(df, var1, var2, lag){

    tidy_cor_test <- ccf(
            df %>% pull(var1),
            df %>% pull(var2),
            lag = lag,
            plot=FALSE
    ) %>%
    broom::tidy() %>%
    mutate(ag_stat = 'do_cross_correlation')
}

tidy_cross_correlation <- function(df, var1, var2, lag, ...){
    df %>%
      group_by_(...) %>%
      do(do_cross_correlation(., var1, var2, lag)) %>%
    mutate(ccf = str_c(var1, ' & ', var2))
}


# exponential decay
do_nls <- function(df, var1, var2){
    tidy_cor_test <- cor.test(
            df %>% pull(var1),
            df %>% pull(var2),
            method = 'pearson'
    ) %>%
    broom::tidy() %>%
    mutate(ag_stat = 'do_cor') %>%
      rename('p_value' = p.value,
             'conf_low' = conf.low,
             'con_high' = conf.high)

    model <- lm(df %>% pull(var2) ~ df %>% pull(var1))

    model %>%
      broom::tidy() %>%
      clean_names() %>%
      mutate(parameter = ifelse(term == '(Intercept)', 'intercept', 'slope')) %>%
      mutate(parameter = str_c('estimate_', parameter)) %>%
      select(parameter, estimate) %>%
      spread(parameter, estimate) %>%
      bind_cols(tidy_cor_test)

}

tidy_nls <- function(df, var1, var2, ...){
    df %>%
      group_by_(...) %>%
      do(do_nls(., var1, var2)) %>%
      quick_sig()
}

# kruskal test (multi group version of ks)
do_kruskal <- function(df, var1, var_group){
    kruskal.test(
        formula(str_c(var1, '~', var_group)),
        df
    ) %>%
    broom::tidy() %>%
    mutate(ag_stat = 'do_kruskal') %>%
    rename('p_value' = p.value)
}

tidy_kruskal <- function(df, var1, var_group, ...){
    df %>%
      group_by(...) %>%
      do(do_kruskal(., var1, var_group)) %>%
      quick_sig()
}

# wilcox test (independent samples, non-parameteric substitute for t test)
do_wilcox_unpaired <- function(df, var1, var2){
    wilcox.test(
            df %>% pull(var1),
            df %>% pull(var2),
            paired = FALSE,
            alternative = 'two.sided'
    ) %>%
    broom::tidy() %>%
    mutate(ag_stat = 'do_wilcox_unpaired')
}

tidy_wilcox_unpaired <- function(df, var1, var2, ...){
    df %>%
      group_by_(...) %>%
      do(do_wilcox_unpaired(., var1, var2))
}

# wilcox test paired (repeated samples, non-parameteric substitute for t test)
do_wilcox_paired <- function(df, var1, var2){
    wilcox.test(
            df %>% pull(var1),
            df %>% pull(var2),
            paired = TRUE,
            alternative = 'two.sided'
    ) %>%
    broom::tidy() %>%
    mutate(ag_stat = 'do_wilcox_unpaired')
}

tidy_wilcox_paired <- function(df, var1, var2, ...){
    df %>%
      group_by_(...) %>%
      do(do_wilcox_paired(., var1, var2))
}


# pairwise wilcox
do_wilcox_pw <- function(df, var1, var_group){
    pairwise.wilcox.test(
        df %>% pull(var1),
        df %>% pull(var_group),
        p.adjust.method="none",
        exact = FALSE
    ) %>%
    broom::tidy() %>%
    mutate(ag_stat = 'do_wilcox_pw') %>%
    rename('p_value' = p.value,
           'group_01' = group1,
           'group_02' = group2)
}

tidy_wilcox_pw <- function(df, var1, var_group, ...){
    df %>%
      group_by(...) %>%
      do(do_wilcox_pw(., var1, var_group)) %>%
      quick_sig()
}


do_sigmoid_pred <- function(df, var_x, var_y){
  y <- df %>% pull(var_y)
  x <- df %>% pull(var_x)

  fit<-NULL
  try(fit<- fit <- nls(y ~ SSlogis(x, Asym, xmid, scal)));

  if(!is.null(fit)){
    x_pred_out <-  data.frame(x = seq(0, 30, length.out = 100))
    y_pred_out <- predict(fit, newdata = x_pred_out)
    return(tibble(x_pred = x_pred_out %>% pull(x), y_pred = as.vector(y_pred_out)))
  } else {
    return(tibble(xpred = NA, y_pred = NA))
  }
}

do_sigmoid_coef <- function(df, var_x, var_y){
  y <- df %>% pull(var_y)
  x <- df %>% pull(var_x)

  fit<-NULL
  try(fit <- nls(y ~ SSlogis(x, Asym, xmid, scal)));

  if(!is.null(fit)){
   coefs <- coef(fit) %>%
       as_tibble() %>%
       mutate(coef = c('asym', 'xmid', 'scal'))

   return(coefs)
  } else {
   coefs <-  tibble(value = c(NA, NA, NA)) %>%
       mutate(coef = c('asym', 'xmid', 'scal'))

   return(coefs)
  }
}

tidy_sigmoid <- function(df, var_x, var_y, ...){
  sigmoid_pred <- df %>%
    group_by_(...) %>%
    do(do_sigmoid_pred(., var_x, var_y))

  sigmoid_coef <- df %>%
    group_by_(...) %>%
    do(do_sigmoid_coef(., var_x, var_y))

  return(list(sigmoid_pred, sigmoid_coef))

}

# anovas
##
## anova one between --------------------------------------------------------------------------------------------------
aov_one_between <-  function(df, dir_output, prefix, var_id, var_dv, var_between){

  anova_model <- aov_ez(
    id = var_id,
    dv = var_dv,
    data = data_plt,
    between = var_between,
    anova_table=list(correction="none", es = "none"))


  anova_summary <-  tibble(
    var_dv = var_dv,
    var_between = var_between,
    df_num = summary(anova_model)[[1]],
    df_den = summary(anova_model)[[2]],
    mse = summary(anova_model)[[3]],
    f = summary(anova_model)[[4]],
    p_value = summary(anova_model)[[5]]
    )

  interaction <- paste("~", var_between, sep = "")

  emm <- emmeans(anova_model, formula(interaction))


  pairs_hsd <- pairs(emm) %>%
    as_tibble() %>%
    separate(contrast, into = c("cont1", "cont2"), "-") %>%
    mutate(cont1 = cont1%>%str_trim(),
           cont2 = cont2%>%str_trim()
           ) %>%
    clean_names()

  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)

  return(list(anova_model, anova_summary, pairs_hsd))

}


# rm anovas
## rm anova one between one within-------------------------------------------------------------------------------------
aov_rm_one_between_one_within <- function(df, dir_output, prefix, var_id, var_dependent, var_between, var_within){
  require(afex)

  anova_model <- aov_ez(
    id = var_id, 
    dv = var_dependent, 
    data = df,
    within = var_within,
    between = var_between,
    anova_table=list(correction="none", es = "none"))
  
  
  label_anova <- c("intercept", var_between, var_within, str_c(var_between, '*', var_within)) %>% as_tibble()
  
  anova_summary <-  summary(anova_model)[[4]][] %>% 
                    as_tibble() %>%
                    bind_cols(label_anova,.) %>%
                    rename("effect" = "value",
                           "ss" = "Sum Sq",
                           "df_num" = "num Df",
                           "ss_error" = "Error SS",
                           "df_den" = "den Df",
                           "f" = "F value",
                           "p_value" = "Pr(>F)") 
  
  anova_summary <- anova_summary %>%
    mutate(stat = 'aov_rm_one_between_one_within',
           var_dv   = var_dependent,
           var_id   = var_id,
           var_between = var_between,
           var_within = var_within
           ) %>%
    select(stat, var_dv, var_id, var_between, var_within, everything()) %>%
    quick_sig()

    interaction <- paste("~", var_within, '*', var_between, sep = "")

    emm <- emmeans(anova_model, formula(interaction))

    pairs_hsd <- pairs(emm) %>%
      as_tibble() %>%
      separate(contrast, into = c('contrast1', 'contrast2'), sep = '-') %>%
      mutate_if(is.character, str_trim) %>%
      separate(contrast1, into = c('contrast1_within', 'contrast1_between'), sep = ' ') %>%
      separate(contrast2, into = c('contrast2_within', 'contrast2_between'), sep = ' ') %>%
      mutate_if(is.character, str_trim) %>%
      rename('p_value' = 'p.value',
             't_ratio' = 't.ratio')  %>%
    quick_sig()

  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)
  
  return(list(anova_model, anova_summary, pairs_hsd))
}

aov_rm_one_between_one_within_grouped <- function(df, dir_output, prefix, var_id, var_dependent, value_dependent, var_between, value_between, var_within, value_within,...){

  anova_summary <- df %>%
    group_by_(..., var_dependent_group = var_dependent, var_between_group = var_between, var_within_group = var_within) %>%
    do(aov_rm_one_between_one_within(., NA, NA, var_id, value_dependent, value_between, value_within)[[2]]) %>%
    select(-var_dv, -var_between, -var_within) %>%
    rename('var_dv' = var_dependent_group, 'var_between' = var_between_group, 'var_within' = var_within_group)

  pairs_hsd <- df %>%
    group_by_(..., var_dependent_group = var_dependent, var_between_group = var_between, var_within_group = var_within) %>%
    do(aov_rm_one_between_one_within(., NA, NA, var_id, value_dependent, value_between, value_within)[[3]]) %>%
    rename('var_dv' = var_dependent_group, 'var_between' = var_between_group, 'var_within' = var_within_group)

  anova_summary <- anova_summary %>%
    mutate(effect = effect %>% str_remove_all('value_'))

  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)

  return(list(anova_summary, pairs_hsd))
}

## rm anova one within ------------------------------------------------------------------------------------------------
aov_rm_one_within <- function(df, dir_output, prefix, var_id, var_dependent, var_within){
  require(afex)

  anova_model <- aov_ez(
    id = var_id, 
    dv = var_dependent, 
    data = df,
    within = var_within,
    anova_table=list(correction="none", es = "none"))
  
  
  label_anova <- c("intercept", var_within) %>% as_tibble()
  
  anova_summary <-  summary(anova_model)[[4]][] %>% 
                    as_tibble() %>%
                    bind_cols(label_anova,.) %>%
                    rename("effect" = "value",
                           "ss" = "Sum Sq",
                           "df_num" = "num Df",
                           "ss_error" = "Error SS",
                           "df_den" = "den Df",
                           "f" = "F value",
                           "p_value" = "Pr(>F)") 
  
  anova_summary <- anova_summary %>%
    mutate(stat = 'aov_rm_one_within',
           var_dv      = var_dependent,
           var_id      = var_id,
           var_within  = var_within) %>%
    select(stat, var_dv, var_id, var_within, everything()) %>%
    quick_sig()
  
  
    interaction <- paste("~", var_within, sep = "")

    emm <- emmeans(anova_model, formula(interaction))

    pairs_hsd <- pairs(emm) %>%
      as_tibble() %>%
      separate(contrast, into = c('contrast1_within', 'contrast2_within'), sep = '-') %>%
      mutate_if(is.character, str_trim) %>%
      rename('p_value' = 'p.value',
             't_ratio' = 't.ratio')  %>%
    quick_sig()
  
  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)

  return(list(anova_model, anova_summary, pairs_hsd))
}

## rm anova one between -----------------------------------------------------------------------------------------------
aov_rm_one_between <- function(df, dir_output, prefix, var_id, var_dependent, var_between){
  require(afex)
  require(emmeans)

  anova_model <- aov_ez(
    id = var_id, 
    dv = var_dependent, 
    data = df,
    between = var_between,
    anova_table=list(correction="none", es = "none"))
  
  anova_summary <- summary(anova_model) %>%
    broom::tidy() %>%
    rename("effect" = "term",
           "df_num" = "num.Df",
           "df_den" = "den.Df",
           "mse" = "MSE",
           "f" = "statistic",
           "p_value" = "p.value"
           )
  anova_summary <- anova_summary %>%
    mutate(stat        = 'aov_rm_one_between',
           var_dv      = var_dependent,
           var_id      = var_id,
           var_between = var_between
    ) %>%
    select(stat, var_dv, var_id, var_between, everything()) %>%
    quick_sig()

    interaction <- paste("~", var_between, sep = "")

    emm <- emmeans(anova_model, formula(interaction))

    pairs_hsd <- pairs(emm) %>%
      as_tibble() %>%
      separate(contrast, into = c('contrast1_between', 'contrast2_between'), sep = '-') %>%
      mutate_if(is.character, str_trim) %>%
      rename('p_value' = 'p.value',
             't_ratio' = 't.ratio')  %>%
    quick_sig()

  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)

  return(list(anova_model, anova_summary, pairs_hsd))
}

aov_rm_one_between_grouped <- function(df, dir_output, prefix, var_id, var_dependent, value_dependent, var_between, value_between,...){

  anova_summary <- df %>%
    group_by_(..., var_dependent_group = var_dependent, var_between_group = var_between) %>%
    do(aov_rm_one_between(., NA, NA, var_id, value_dependent, value_between)[[2]]) %>%
    select(-var_dv, -var_between) %>%
    rename('var_dv' = var_dependent_group, 'var_between' = var_between_group)

  pairs_hsd <- df %>%
    group_by_(..., var_dependent_group = var_dependent, var_between_group = var_between) %>%
    do(aov_rm_one_between(., NA, NA, var_id, value_dependent, value_between)[[3]]) %>%
    rename('var_dv' = var_dependent_group, 'var_between' = var_between_group)

  anova_summary <- anova_summary %>%
    mutate(effect = effect %>% str_remove_all('value_'))

  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)

  return(list(anova_summary, pairs_hsd))
}

## rm anova two within ------------------------------------------------------------------------------------------------
aov_rm_two_within <- function(df, dir_output, prefix, var_id, var_dependent, var_within1, var_within2){
  require(afex)

  anova_model <- aov_ez(
    id = var_id, 
    dv = var_dependent, 
    data = df,
    within = c(var_within1, var_within2),
    anova_table=list(correction="none", es = "none"))
  
    label_anova <- c("intercept",
                   var_within1,
                   var_within2,
                   str_c(var_within1, '*', var_within2)
                   ) %>%
    as_tibble()

    anova_summary <-  summary(anova_model)[[4]][] %>%
                    as_tibble() %>%
                    bind_cols(label_anova,.) %>%
                    rename("effect" = "value",
                           "ss" = "Sum Sq",
                           "df_num" = "num Df",
                           "ss_error" = "Error SS",
                           "df_den" = "den Df",
                           "f" = "F value",
                           "p_value" = "Pr(>F)")

  
  anova_summary <- anova_summary %>%
    mutate(stat = 'aov_rm_two_within',
           var_dv   = var_dependent,
           var_id   = var_id,
           var_within1 = var_within1,
           var_within2 = var_within2
    ) %>%
    select(stat, var_dv, everything()) %>%
    quick_sig()

    interaction <- paste("~", var_within1, '*', var_within2, sep = "")

    emm <- emmeans(anova_model, formula(interaction))

    pairs_hsd <- pairs(emm) %>%
      as_tibble() %>%
      separate(contrast, into = c('contrast1', 'contrast2'), sep = '-') %>%
      mutate_if(is.character, str_trim) %>%
      separate(contrast1, into = c('contrast1_within1', 'contrast1_within2'), sep = ' ') %>%
      separate(contrast2, into = c('contrast2_within1', 'contrast2_within2'), sep = ' ') %>%
      mutate_if(is.character, str_trim) %>%
      rename('p_value' = 'p.value',
             't_ratio' = 't.ratio')  %>%
    quick_sig()

  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)

  return(list(anova_model, anova_summary, pairs_hsd))

}

## rm anova two within ------------------------------------------------------------------------------------------------
aov_rm_three_within <- function(df, dir_output, prefix, var_id, var_dependent, var_within1, var_within2, var_within3){
  require(afex)

  anova_model <- aov_ez(
    id = var_id,
    dv = var_dependent,
    data = df,
    within = c(var_within1, var_within2, var_within3),
    anova_table=list(correction="none", es = "none"))

    label_anova <- c("intercept",
                   var_within1,
                   var_within2,
                   var_within3,
                   str_c(var_within1, '*', var_within2),
                   str_c(var_within1, '*', var_within3),
                   str_c(var_within2, '*', var_within3),
                   str_c(var_within1, '*', var_within2, '*', var_within3)
                   ) %>%
    as_tibble()

    anova_summary <-  summary(anova_model)[[4]][] %>%
                    as_tibble() %>%
                    bind_cols(label_anova,.) %>%
                    rename("effect" = "value",
                           "ss" = "Sum Sq",
                           "df_num" = "num Df",
                           "ss_error" = "Error SS",
                           "df_den" = "den Df",
                           "f" = "F value",
                           "p_value" = "Pr(>F)")


  anova_summary <- anova_summary %>%
    mutate(stat = 'aov_rm_three_within',
           var_dv   = var_dependent,
           var_id   = var_id,
           var_within1 = var_within1,
           var_within2 = var_within2
    ) %>%
    select(stat, var_dv, everything()) %>%
    quick_sig()

    interaction <- paste("~", var_within1, '*', var_within2, '*', var_within3, sep = "")

    emm <- emmeans(anova_model, formula(interaction))

    pairs_hsd <- pairs(emm) %>%
      as_tibble() %>%
      separate(contrast, into = c('contrast1', 'contrast2'), sep = '-') %>%
      mutate_if(is.character, str_trim) %>%
      separate(contrast1, into = c('contrast1_within1', 'contrast1_within2', 'contrast1_within3'), sep = ' ') %>%
      separate(contrast2, into = c('contrast2_within1', 'contrast2_within2', 'contrast2_within3'), sep = ' ') %>%
      mutate_if(is.character, str_trim) %>%
      rename('p_value' = 'p.value',
             't_ratio' = 't.ratio')  %>%
    quick_sig()

  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)

  return(list(anova_model, anova_summary, pairs_hsd))

}

## rm anova one between two within -------------------------------------------------------------------------------------
aov_rm_one_between_two_within <- function(df, dir_output, prefix, var_id, var_dependent, var_between, var_within1, var_within2){
  require(afex)
  require(emmeans)

  anova_model <- aov_ez(
    id = var_id,
    dv = var_dependent,
    data = df,
    between = var_between,
    within = c(var_within1, var_within2),
    anova_table=list(correction="none", es = "none"))

  label_anova <- c("intercept",
                   var_between,
                   var_within1,
                   str_c(var_between, '*', var_within1),
                   var_within2,
                   str_c(var_between, '*', var_within2),
                   str_c(var_within1, '*', var_within2),
                   str_c(var_between, '*', var_within1, '*', var_within2)
                   ) %>%
    as_tibble()

  anova_summary <-  summary(anova_model)[[4]][] %>%
                    as_tibble() %>%
                    bind_cols(label_anova,.) %>%
                    rename("effect" = "value",
                           "ss" = "Sum Sq",
                           "df_num" = "num Df",
                           "ss_error" = "Error SS",
                           "df_den" = "den Df",
                           "f" = "F value",
                           "p_value" = "Pr(>F)")

  anova_summary <- anova_summary %>%
    mutate(stat = 'aov_rm_one_between_two_within',
           var_dv   = var_dependent,
           var_id   = var_id,
           var_between = var_between,
           var_within1 = var_within1,
           var_within2 = var_within2) %>%
    select(stat, var_dv, var_id, var_between, var_within1, var_within2, everything())  %>%
    quick_sig()

  interaction <- paste("~", var_between, '*', var_within1, '*', var_within2, sep = "")

  emm <- emmeans(anova_model, formula(interaction))

  pairs_hsd <- pairs(emm) %>%
    as_tibble() %>%
    separate(contrast, into = c('contrast1', 'contrast2'), sep = '-') %>%
    mutate_if(is.character, str_trim) %>%
    separate(contrast1, into = c('contrast1_between', 'contrast1_within1', 'contrast1_within2'), sep = ' ') %>%
    separate(contrast2, into = c('contrast2_between', 'contrast2_within1', 'contrast2_within2'), sep = ' ') %>%
    mutate_if(is.character, str_trim) %>%
    rename('p_value' = 'p.value',
           't_ratio' = 't.ratio')  %>%
    quick_sig()

  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)

  return(list(anova_model, anova_summary, pairs_hsd))
}
