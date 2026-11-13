#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("datetime", quietly = TRUE)) {
  install.packages("datetime")
}

library(datetime)
library(lubridate)

exe <- "./stress_social.exe"
file_name_prefix <- "sim_stress_social"

npatches_vals = c(5000) # number of patches
n_vals = c(1) # num individuals per patch
max_time_vals = c(1000000) # max time
risk <- c(0.05, 0.10)
autocorrelation <- c(0.0, 0.3, 0.9)
md_vals = c(1.0) # weight of damage-related mortality
mv_vals = c(0.0) # weight of vigilance-related mortality
p_mig_vals = c(1.0) # migration probability
p_attack_vals = c(0.5) # probability of being attacked when predator present
init_v_vals = c(0) # Initial vigilance
init_stress_vals = c(0.4) # Initial stress hormone level
g_vals = c(1.0) # Damage removal per timestep
k_vals = c(1.0) # Increase in damage due to hormone != optimum

ctr <- 0

fecundity_power <- 1.0
hmax <- 10.0
vigilance_on <- 0
run_end_assay <- 0
assay_n_individuals <- 0
assay_pre_time <- 0 
assay_post_time <- 0 

# date time
current_time <- now()
formatted_time <- format(current_time, "%Y%m%d_%H%M%S")

for (npatches_i in npatches_vals) {
    for (n_i in n_vals) {
        for (max_time_i in max_time_vals) {
            for (risk_i in risk) {
                for (autocorrelation_i in autocorrelation) {
                    s_np_i <- risk_i * (1 - autocorrelation_i)
                    s_p_i  <- (1 - risk_i) * (1 - autocorrelation_i)
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
                                                    output_file_name <- paste0(
                                                            file_name_prefix, "_", 
                                                            ctr, "_", 
                                                            formatted_time)

                                                    # Write the command line with parameters
                                                    writeLines(text = paste0(
                                                        exe," ",
                                                        output_file_name," ",
                                                        npatches_i," ",
                                                        n_i," ",
                                                        max_time_i," ",
                                                        s_np_i," ",
                                                        s_p_i, " ",
                                                        md_i," ",
                                                        mv_i," ", 
                                                        p_mig_i," ",
                                                        p_attack_i, " ",
                                                        fecundity_power," ",
                                                        hmax," ",
                                                        init_v_i," ",
                                                        init_stress_i, " ",
                                                        g_i, " ",
                                                        k_i, " ",
                                                        vigilance_on," ",
                                                        run_end_assay," ",
                                                        assay_n_individuals," ",
                                                        assay_pre_time," ",
                                                        assay_post_time," "
                                                        ))
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

