#!/usr/bin/env Rscript

exe <- "./stress_social.exe"
file_name_prefix <- "sim_stress_social"

s_np = c(0.1,0.5,0.9)
s_p = c(0.1,0.5,0.9)

npatches = 200
n = 20

ctr <- 0



# TODO
# make this date time dependent using something like datetime

for (s_np_i in s_np) 
{
    for (s_p_i in s_p) 
    {
        ctr <- ctr + 1
        writeLines(text = paste0(exe," ",file_name_prefix
                        ,"_",ctr," ",s_np_i
                        , " ", s_p_i, " ",npatches, " ",n))
    }
}


