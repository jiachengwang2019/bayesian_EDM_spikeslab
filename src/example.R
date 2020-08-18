#  #######################################################################
#       File-Name:      example.R
#       Version:        R 3.6.2
#       Date:           Aug, 18 2020
#       Author:         Jiacheng Wang <Jiacheng.Wang@nbcuni.com>
#       Purpose:        Example to use find_champion_bayesian_ss function
#       Input Files:    NONE
#       Output Files:   NONE
#       Data Output:    NONE
#       Previous files: NONE
#       Dependencies:   NONE
#       Required by:    NONE
#       Status:         IN PROGRESS
#       Machine:        NBCU laptop
#  #######################################################################
library(here)
source(here::here("src","find_champion_bayesian_ss.R"))
source(here::here("src","parameters.R"))

options(warn = -1,digits = 3,verbose = F,error = NULL)
raw_data_AFL = create_raw_data(inputFileName = inFileName,
                               inputFilePath = inFilePath_AFL,
                               keep_cols,
                               group_by_cols)

case_data_SalesPrime <- create_aggregated_data(raw_data_AFL, 
                                               filter_text = "HH_NT_Daypart == 'SalesPrime' & is.na(SC3_Impressions)==F",
                                               network = "Cable",
                                               agg_timescale = "Week",
                                               agg_group_by_cols = weekly_group_by_cols_new)

# Setup the weekly/daily regressors
wkly_regressors = c(
  "Jan","Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Oct", "Nov", "Dec",
  "Easter_week_ind", "Memorial_Day_week_ind", "Independence_Day_week_ind", 
  "Halloween_week_ind", "Thanksgiving_week_ind", "Thanksgiving_week_aft_ind", 
  "Christmas_week_bfr_ind","Christmas_week_ind","New_Years_Eve_week_ind"
)
regressors_intercept_weekly = c(wkly_regressors, "trend", "intercept")
daily_regressors = c(
  "Jan","Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Oct", "Nov", "Dec",
  "Sun","Mon","Tue","Thu","Fri","Sat",
  "Easter_ind", "Memorial_Day_ind", "Independence_Day_ind", "Halloween_ind", 
  "Thanksgiving_ind", "Christmas_ind","New_Years_Eve_ind"
)
regressors_intercept_daily = c(daily_regressors, "trend", "intercept")

#:::::::::::::::::::::::::::::::::::::::::::: Test find_champion_bayesian_ss function
res <- find_champion_bayesian_ss(data = case_data_SalesPrime,
                                 agg_timescale = 'Week',
                                 OOS_start = as.Date('2018-01-01'),
                                 regressors = regressors_intercept_weekly,
                                 nr_samples = 50,
                                 nr_burnin = 0)

res$champion_model
res$champion_result
res$all_model_details

