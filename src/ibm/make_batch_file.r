#!/usr/bin/env Rscript

library(datetime)
library(lubridate)

exe <- "./stress_social.exe"
file_name_prefix <- "sim_stress_social"

npatches_vals = c(200) # number of patches
n_vals = c(20) # num individuals per patch
max_time_vals = c(30000) # max time
s_np = c(0.1) # switch rate NP to P
s_p = c(0.1) # switch rate P to NP 
md_vals = c(1.0) # weight of damage-related mortality
mv_vals = c(1.0, 5.0, 10.0) # weight of vigilance-related mortality
p_mig_vals = c(0.1) # migration probability
p_attack_vals = c(0.1) # probability of being attacked when predator present
init_v_vals = c(0.3) # Initial vigilance
init_stress_vals = c(0.4) # Initial stress hormone level
g_vals = c(0.1) # Damage removal per timestep
k_vals = c(0.1) # Increase in damage due to hormone != optimum

ctr <- 0

# date time
current_time <- now()
formatted_time <- format(current_time, "%Y%m%d_%H%M%S")

for (npatches_i in npatches_vals) {
    for (n_i in n_vals) {
        for (max_time_i in max_time_vals) {
            for (s_np_i in s_np) {
                for (s_p_i in s_p) {
                    for (md_i in md_vals) {
                        for (mv_i in mv_vals) {
                            for (p_mig_i in p_mig_vals) {
                                for (p_attack_i in p_attack_vals) {
                                    for (init_v_i in init_v_vals) {
                                        for (init_stress_i in init_stress_vals) {
                                            for (g_i in g_vals) {
                                                for (k_i in k_vals) {
                                                    
                                                    ctr <- ctr + 1

                                                    # Output file name
                                                    output_file_name <- paste0(file_name_prefix, "_", ctr, "_", formatted_time)

                                                    # Write the command line with parameters
                                                    writeLines(text = paste0(
                                                        exe, " ", output_file_name, " ", npatches_i, 
                                                        " ", n_i, " ", max_time_i, " ", s_np_i, 
                                                        " ", s_p_i, " ", md_i, " ", mv_i, 
                                                        " ", p_mig_i, " ", p_attack_i, " ", init_v_i, 
                                                        " ", init_stress_i, " ", g_i, " ", k_i))
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

