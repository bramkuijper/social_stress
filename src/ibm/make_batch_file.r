#!/usr/bin/env Rscript

library(datetime)
library(lubridate)

exe <- "./stress_social.exe"
file_name_prefix <- "sim_stress_social"

s_np = c(0.1,0.5,0.9)
s_p = c(0.1,0.5,0.9)
mv_vals = c(0.1, 0.5, 1.0) #new values for mv

npatches = 200
n = 20

ctr <- 0

# date time addition

current_time <- now()
formatted_time <- format(current_time, "%Y%m%d_%H%M%S")

for (s_np_i in s_np) 
{
    for (s_p_i in s_p) {
        for (mv_i in mv_vals) { 

        
        ctr <- ctr + 1

        # output file name
        output_file_name <- paste0(file_name_prefix,"_",ctr,"_",formatted_time) 
        # assuming it needs a file type? There wasn't one before

        writeLines(text = paste0(
                        exe," ",output_file_name," ",s_np_i,
                         " ", s_p_i, " ",npatches, " ",n, " ", mv_i))
        }
    }
}


