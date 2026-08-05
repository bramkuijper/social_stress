#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

enum Sex
{
    female = 0,
    male = 1
};

enum PatchType
{
    NP = 0,
    P = 1
};

// all parameters with their default values
class Parameters
{
    public:

        unsigned data_output_interval{1};

        // population size (no sex differences)
        unsigned n{20};

        // dispersal
        double p_mig{0.1};

        unsigned npatches{100}; 

        unsigned max_time{30000};

        // switch rates between predator present in patch
        // vs absent
        double s[2]{0.5,0.5};

        // attack probability in a patch where a predator
        // is present

        double p_attack{0.1};

        // power of how fecundity decreases with vigilance
        double fecundity_power{1.0};
        
        // power controlling how damage reduces fecundity.
        // Added for Taborsky-validation runs so that the cost of elevated
        // stress hormone acts through reproduction rather than mortality,
        // matching the fitness trade-off in the Taborsky stress model.
        double damage_fecundity_power{1.5};

        // min max hormone
        double hmin{0.0};
        double hmax{10.0};

        // Maximum damage level.
        // Set to 1.0 for Taborsky validation so that damage uses the same
        // [0,1] scale as hormone
        double dmax{1.0};
        
        // Power controlling the protective effect of hormone during an attack.
        // Set to 1.0 for Taborsky validation
        double survival_power{1.0};

        double init_v{0.0};
        double init_stress_hormone_level{0.0};
        double init_removal{0.1};
        
        // TABORSKY VALIDATION:
        // Maximum stress-induced hormone influx per timestep during
        // an active post-stressor response.
        double stress_influx_max{0.25};
        
        // Number of timesteps for which stress-induced hormone production
        // can continue after an attack. Taborsky Box 3 used 75.
        unsigned tmax_stress_influx{75};
        
        // Initial value of h1_S, which controls negative feedback on
        // stress-induced hormone production.
        double init_h1_S{0.0};

        // base name for the file
        std::string file_name{"sim_stress_social"};

        // TABORSKY VALIDATION:
        // mutation rates aligned with Taborsky stress model
        double mu_baseline{0.0005};
        double mu_stress_influx{0.0};
        double mu_h1_S{0.005};
        double mu_vigilance_influx{0.0}; // Mutation rate for vigilance-driven stress influx - set to 0.0 in no-vigilance benchmark runs
        double mu_removal{0.0005};
        double mu_v{0.0}; // Mutation rate for baseline vigilance - set to 0.0 in no-vigilance benchmark runs so vigilance can't evolve
        // Standard deviation of mutational effects
        // Taborsky validation uses 0.10
        double sdmu{0.10};

        // mortality rates 
        double m0{0.001}; // 1/1000 mortality
        double md{1.0}; // weighting of damage-related mortality
        double mv{1.0}; // weighting of vigilance-investment-related mortality
        
        // damage-related things
        double g{0.1}; // removal of damage per timestep
        double k{0.1}; // increase in damage due to hormone != optimum
        // Hormone optimum used by the original social_stress damage functoin
        // Not used in the Taborsky-validation damage equation but retaining in case it links with McN vigilance links
        double theta_hormone{1}; // optimal hormone level
        
        // TABORSKY VALIDATION:
        // Initial stress-independent hormone influx.
        double init_baseline_influx{0.05};
        
        // Starting hormone level is an evolvable trait in the
        // Taborsky stress model.
        double init_hstart{0.5};
        
        // Mutation probability for hstart.
        double mu_hstart{0.0005};
        
        // Minimum clearance allowed in the Taborsky model.
        double min_removal{0.01};
        
        // whether vigilance is enabled in this simulation - existing behaviour is default
        bool vigilance{true};
        
        // whether to run the controlled end-of-simulation hormone assay
        bool run_end_assay{false};

        // number of evolved individuals sampled for the assay
        unsigned assay_n_individuals{10};
        
        // number of timesteps before the standardised attack
        unsigned assay_pre_time{25};
        
        // number of timesteps after the standardised attack
        unsigned assay_post_time{100};
};

#endif
