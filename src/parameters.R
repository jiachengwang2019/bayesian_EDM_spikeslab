#  #######################################################################
#       File-Name:      parameters.R
#       Version:        R 3.5.1
#       Date:           Jan 09, 2019
#       Author:         Soudeep Deb <Soudeep.Deb@nbcuni.com>
#       Purpose:        defines parameters that will be used in the pacing models
#       Input Files:    NONE
#       Output Files:   NONE
#       Data Output:    NONE
#       Previous files: NONE
#       Dependencies:   NONE
#       Required by:    MAIN.R
#       Status:         APPROVED
#       Machine:        NBCU laptop
#  #######################################################################
# :::::::::DATA-RELATED PARAMETERS ::::::::::::::::::::::::::::::::::::::::::::
# ::::::: GENERAL PATH & DATA
inFilePath_AFL <- "/mnt/nbcu-ds-linear-forecasting/processed/pacing/19W019/USA_AFL_0325_190326_19W019.csv"
inFileName <- "USA_AFL_0325_190326_19W019"
# :::::::: GROUP_BY VARIABLES
group_by_cols <- c("Source","Network", "Demo", "HH_NT_Daypart","Broadcast_Year","Cable_Qtr","Week","Date",
                   "Half_Hr","Start_Time","Show_Name","program_type","FirstAiring",
                   "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug","Sep","Oct", "Nov", "Dec", 
                   "New_Years_Eve_ind","Memorial_Day_week_ind","Memorial_Day_ind","Halloween_ind",              
                   "Thanksgiving_ind","Christmas_week_ind","Christmas_ind","Christmas_week_bfr_ind",    
                   "Halloween_week_ind","Thanksgiving_week_ind","Independence_Day_week_ind","Easter_week_ind",             
                   "New_Years_Eve_week_ind","Easter_ind","Independence_Day_ind","Thanksgiving_week_aft_ind","Season_num","Genre","Episode")
hourly_group_by_cols <- c("Source","Network", "Demo", "HH_NT_Daypart","Broadcast_Year",
                          "Cable_Qtr","Week","Date","FirstAiring","Show_Name","program_type",
                          "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                          "Jul", "Aug","Sep","Oct", "Nov", "Dec", 
                          "New_Years_Eve_ind","Memorial_Day_week_ind","Memorial_Day_ind","Halloween_ind",              
                          "Thanksgiving_ind","Christmas_week_ind","Christmas_ind","Christmas_week_bfr_ind",    
                          "Halloween_week_ind","Thanksgiving_week_ind","Independence_Day_week_ind","Easter_week_ind",             
                          "New_Years_Eve_week_ind","Easter_ind","Independence_Day_ind","Thanksgiving_week_aft_ind","Season_num","Genre","Episode")
daily_group_by_cols <- c("Source","Network","Demo","HH_NT_Daypart","Broadcast_Year",
                         "Cable_Qtr","Week","Date","Jan","Feb","Mar","Apr","May","Jun",
                         "Jul","Aug","Sep","Oct","Nov","Dec", 
                         "New_Years_Eve_ind","Memorial_Day_week_ind","Memorial_Day_ind","Halloween_ind",              
                         "Thanksgiving_ind","Christmas_week_ind","Christmas_ind","Christmas_week_bfr_ind",    
                         "Halloween_week_ind","Thanksgiving_week_ind","Independence_Day_week_ind","Easter_week_ind",             
                         "New_Years_Eve_week_ind","Easter_ind","Independence_Day_ind","Thanksgiving_week_aft_ind","Season_num","Genre","Episode")
daily_group_by_cols_new <- c("Source","Network","Demo","HH_NT_Daypart","Broadcast_Year",
                             "Cable_Qtr","Week","Date","Jan","Feb","Mar","Apr","May","Jun",
                             "Jul","Aug","Sep","Oct","Nov","Dec", 
                             "New_Years_Eve_ind","Memorial_Day_week_ind","Memorial_Day_ind","Halloween_ind",              
                             "Thanksgiving_ind","Christmas_week_ind","Christmas_ind","Christmas_week_bfr_ind",    
                             "Halloween_week_ind","Thanksgiving_week_ind","Independence_Day_week_ind","Easter_week_ind",             
                             "New_Years_Eve_week_ind","Easter_ind","Independence_Day_ind","Thanksgiving_week_aft_ind")
weekly_group_by_cols <- c("Source","Network","Demo","HH_NT_Daypart","Broadcast_Year",
                          "Cable_Qtr","Week","Jan","Feb","Mar","Apr","May","Jun",
                          "Jul","Aug","Sep","Oct","Nov","Dec", 
                          "Memorial_Day_week_ind","Christmas_week_ind","Christmas_week_bfr_ind",    
                          "Halloween_week_ind","Thanksgiving_week_ind","Independence_Day_week_ind","Easter_week_ind",             
                          "New_Years_Eve_week_ind","Thanksgiving_week_aft_ind","Season_num","Genre","Episode")
weekly_group_by_cols_new <- c("Source","Network","Demo","HH_NT_Daypart","Broadcast_Year",
                              "Cable_Qtr","Week","Jan","Feb","Mar","Apr","May","Jun",
                              "Jul","Aug","Sep","Oct","Nov","Dec", 
                              "Memorial_Day_week_ind","Christmas_week_ind","Christmas_week_bfr_ind",    
                              "Halloween_week_ind","Thanksgiving_week_ind","Independence_Day_week_ind","Easter_week_ind",             
                              "New_Years_Eve_week_ind","Thanksgiving_week_aft_ind")
# :::::::::: VARIABLES TO KEEP IN THE INITIAL DATA FILTERING
keep_cols <- c("Source","Network", "Demo", "Broadcast_Year", "Cable_Qtr", "Date", "Week",
               "Show_Name", "NBC_Show_Name", "ShowName_Orig","program_type",
               "Start_Time", "End_Time", "Half_Hr", "HH_NT_Daypart","FirstAiring",
               "Tot_UE", "LS_Imps", "Nielsen_LS_Rating", "LS_Dur",
               "SC3_Impressions", "Nielsen_SC3_Rating", "SC3_C_Dur",
               "Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
               "Sun","Mon","Tue","Wed","Thu","Fri","Sat", 
               "New_Years_Eve_ind","Memorial_Day_week_ind","Memorial_Day_ind","Halloween_ind",              
               "Thanksgiving_ind","Christmas_week_ind","Christmas_ind","Christmas_week_bfr_ind",    
               "Halloween_week_ind","Thanksgiving_week_ind","Independence_Day_week_ind","Easter_week_ind",             
               "New_Years_Eve_week_ind","Easter_ind","Independence_Day_ind","Thanksgiving_week_aft_ind","Season_num","Genre","Episode")
# ::::::::: VARIABLES TO KEEP IN THE OUTPUT DATA
output_cols <- c("Network", "Demo", "HH_NT_Daypart", "Broadcast_Year", "Cable_Qtr", "Week", "Show_Name",
                 "SC3_Impressions", "C3_Rating","SC3_C_Dur","LS_Rating","LS_Dur","Stream","Predict")

