// heart of the stress social code

#include <cassert>
#include <iostream>
#include <cmath>
#include <numeric>
#include <fstream>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <utility>

#include "stress_social.hpp"
#include "patch.hpp"
#include "individual.hpp"

// vigilance can be switched on or off depending on input for sim run
inline double effective_vigilance(
          Individual const & ind,
          bool const vigilance_enabled)
{
    if (!vigilance_enabled)
    {
      return 0.0;
    }
    
    double base_v = 0.5 * (ind.v[0] + ind.v[1]);

    if (base_v < 0.0) { base_v = 0.0; }
    if (base_v > 1.0) { base_v = 1.0; }

    return base_v;
}


// constructor function
StressSocial::StressSocial(Parameters const &parvals) :
    param{parvals}  // Random number generator initialisation 
    ,rd{} // initialize random device
    ,seed{rd()} // initialize seed
    ,rng_r{seed} // initialize the random number generator
    ,uniform{0.0,1.0} // initialize the uniform dist between 0 and 1
    ,patch_sampler{0, param.npatches - 1} // initialize uniform distribution to sample patch indices from
    ,take_random_breeder{0, param.n - 1} // initialize uniform distribution to sample patch indices from
    ,metapopulation(param.npatches, Patch(param)) // initialize the metapopulation
    ,data_file{param.file_name} // File where output is written
    ,last_total_global_fecundity{0.0} // total fecundity across all patches last timestep
{
    write_data_headers();

    // make some patches P and some NP
    initialize_patches();

    // now run the thing (with EG FIXES to update predator presence/absence for each patch before
    // writing of distribution to new file at time_step == 0 
    for (time_step = 0; time_step < param.max_time; ++time_step)
    {
          // at start of first time step, write out all individuals - comment out if other section for first "iffy" timestep is included
          // if (time_step == 0) {
          //    write_distribution();
          //    }
          
          // reset counters at the start of each timestep
            n_attacked = 0; 
            n_death_damage = 0;
            n_death_predator = 0;
            sum_damage_at_damage_death = 0.0;

          // TABORSKY VALIDATION:
          // Match the order of events in the Taborsky stress model:
          //
          // 1. environment changes
          // 2. baseline hormone dynamics
          // 3. predator attack
          // 4. stress-induced hormone response
          // 5. background mortality
          // 6. damage update
          // 7. fecundity/reproduction
          switch_predator_status();
          
          update_baseline_hormone();
          
          predator_visit();
          
          update_stress_response();
          
          survive_damage_vigilance();
          
          update_damage();
          
          reproduce();
        
           // error checking: ntotal should always be >= each death count - simplified from previous version
            assert(param.n * param.npatches >= n_death_damage); // total pop >= deaths from damage
            assert(param.n * param.npatches >= n_death_predator); // total pop >= deaths from predation
            assert(param.n * param.npatches >= n_death_damage + n_death_predator); // total pop >= total deaths
        
            if (time_step % param.data_output_interval == 0)
            {
                write_data();
                // write_distribution(); // debug only
            }
    } // end evolutionary simulation
    
    // Write every individual in the final evolved population
    // Done before optional assay so that the file represents the population 
    // immediately after evolution
    
    write_final_individuals();
    
    // Optional controlled hormone assay using evolved individuals
    if (param.run_end_assay)
    {
        run_end_hormone_assay();
    }

    // Write simulation parameters to the main output file
    write_parameters();
    
} // end StressSocial constructor

// Write the complete diploid genotype of every individual in the population after the evolutionary simulation has finished
// Function called once so only contains the final evolved generation
void StressSocial::write_final_individuals()
{
    namespace fs = std::filesystem;

    // Split the main output path into its directory and filename.
    fs::path main_output_path{param.file_name};

    fs::path parent_directory =
        main_output_path.parent_path();

    // If param.file_name contains no directory component,
    // create the individuals folder in the current directory.
    if (parent_directory.empty())
    {
        parent_directory = ".";
    }

    // Create a subfolder called "individuals" alongside
    // the normal simulation output files.
    fs::path individuals_directory =
        parent_directory / "individuals";

    std::error_code directory_error;

    fs::create_directories(
        individuals_directory,
        directory_error
    );

    if (directory_error)
    {
        std::cerr
            << "Could not create individuals output directory: "
            << individuals_directory
            << "\nReason: "
            << directory_error.message()
            << "\n";

        return;
    }

    // Use the normal simulation filename, with "_individuals"
    // added to the end.
    fs::path individuals_path =
        individuals_directory /
        (main_output_path.filename().string() + "_individuals");

    individuals_file.open(individuals_path);

    if (!individuals_file)
    {
        std::cerr
            << "Could not open final-individual output file: "
            << individuals_path
            << "\n";

        return;
    }

    // Each row represents one individual in the final population.
    // The two alleles at each diploid locus are written separately.
    individuals_file
        << "final_time_step;"
        << "patch_index;"
        << "breeder_index;"
        << "v0;"
        << "v1;"
        << "baseline_influx0;"
        << "baseline_influx1;"
        << "stress_influx0;"
        << "stress_influx1;"
        << "h1_S0;"
        << "h1_S1;"
        << "hstart0;"
        << "hstart1;"
        << "vigilance_influx0;"
        << "vigilance_influx1;"
        << "removal0;"
        << "removal1;"
        << "stress_hormone;"
        << "damage"
        << '\n';

    // Loop once over every individual in the final population.
    for (unsigned patch_idx = 0;
         patch_idx < metapopulation.size();
         ++patch_idx)
    {
        const Patch &patch =
            metapopulation[patch_idx];

        for (unsigned breeder_idx = 0;
             breeder_idx < patch.breeders.size();
             ++breeder_idx)
        {
            const Individual &individual =
                patch.breeders[breeder_idx];

            individuals_file
                << (param.max_time > 0
                        ? param.max_time - 1
                        : 0)
                << ";"
                << patch_idx << ";"
                << breeder_idx << ";"
                << individual.v[0] << ";"
                << individual.v[1] << ";"
                << individual.baseline_influx[0] << ";"
                << individual.baseline_influx[1] << ";"
                << individual.stress_influx[0] << ";"
                << individual.stress_influx[1] << ";"
                << individual.h1_S[0] << ";"
                << individual.h1_S[1] << ";"
                << individual.hstart[0] << ";"
                << individual.hstart[1] << ";"
                << individual.vigilance_influx[0] << ";"
                << individual.vigilance_influx[1] << ";"
                << individual.removal[0] << ";"
                << individual.removal[1] << ";"
                << individual.stress_hormone << ";"
                << individual.damage
                << '\n';
        }
    }

    // The final population is written only once, so close
    // the file immediately after completing the output.
    individuals_file.close();
}

// go over all the patches and initialize them as type NP or P
void StressSocial::initialize_patches()
{
    // calculate probability of encountering a predator on a patch using switch rates
    double prob_P = param.s[NP] / (param.s[NP] + param.s[P]);

    for (auto patch_iterator = metapopulation.begin();
            patch_iterator != metapopulation.end();
            ++patch_iterator)
    {
        patch_iterator->predator_patch = uniform(rng_r) < prob_P; // If draw number lower than prob_P, then P = TRUE
    }
}

// update whether each patch has a predator, once per timestep
// Uses s[NP] (NP -> P) and s[P] (P -> NP)

void StressSocial::switch_predator_status()
{

    for (auto metapop_iter = metapopulation.begin();
         metapop_iter != metapopulation.end();
         ++metapop_iter)
    {
        if (metapop_iter->predator_patch == false)
        {
            // currently NP; can switch to P with probability s[NP]
            if (uniform(rng_r) < param.s[NP])
            {
                metapop_iter->predator_patch = true;
            }
        }
        else
        {
            // currently P; can switch to NP with probability s[P]
            if (uniform(rng_r) < param.s[P])
            {
                metapop_iter->predator_patch = false;
            }
        }
    }
}

// print out the distribution/values of all the individuals
// creation of separate debug file
void StressSocial::write_distribution()
{
    // Separate debug file: same base name +"_distribution"
    std::ofstream data_file2{param.file_name + "_distribution"};
    
    if (!data_file2) {
        std::cerr << "Could not open distribution output file: "
                  << param.file_name + "_distribution" << "\n";
        return;
    }

    // Header row
    data_file2 << "time_step;"
               << "patch_index;"
               << "breeder_index;"
               << "is_alive;"
               << "is_attacked;"
               << "predator_patch;"
               << "V;"
               << "v0;v1;"
               << "baseline_influx0;baseline_influx1;"
               << "stress_influx0;stress_influx1;"
               << "h1_S0;h1_S1;"
               << "hstart0;hstart1;"
               << "vigilance_influx0;vigilance_influx1;"
               << "removal0;removal1;"
               << "damage;"
               << "stress_hormone"
               << '\n';

    // Loop over all patches and all breeders
    for (unsigned patch_idx = 0; patch_idx < metapopulation.size(); ++patch_idx) {
        const Patch &patch = metapopulation[patch_idx];

        for (unsigned breeder_idx = 0;
             breeder_idx < patch.breeders.size();
             ++breeder_idx) {

            const Individual &ind = patch.breeders[breeder_idx];

            data_file2 << time_step << ";"
                       << patch_idx << ";"
                       << breeder_idx << ";"
                       << ind.is_alive << ";"
                       << ind.is_attacked << ";"
                       << patch.predator_patch << ";"
                       << patch.V << ";"
                       << ind.v[0] << ";"
                       << ind.v[1] << ";"
                       << ind.baseline_influx[0] << ";"
                       << ind.baseline_influx[1] << ";"
                       << ind.stress_influx[0] << ";"
                       << ind.stress_influx[1] << ";"
                       << ind.h1_S[0] << ";"
                       << ind.h1_S[1] << ";"
                       << ind.hstart[0] << ";"
                       << ind.hstart[1] << ";"
                       << ind.vigilance_influx[0] << ";"
                       << ind.vigilance_influx[1] << ";"
                       << ind.removal[0] << ";"
                       << ind.removal[1] << ";"
                       << ind.damage << ";"
                       << ind.stress_hormone
                       << '\n';
        }
    }

}

// Run a controlled hormone-response assay on a random sample of
// evolved individuals after the main simulation has finished.
void StressSocial::run_end_hormone_assay()
{
    // Store every possible patch/breeder location.
    std::vector<std::pair<unsigned, unsigned>> individual_locations;

    for (unsigned patch_idx = 0;
         patch_idx < metapopulation.size();
         ++patch_idx)
    {
        for (unsigned breeder_idx = 0;
             breeder_idx < metapopulation[patch_idx].breeders.size();
             ++breeder_idx)
        {
            individual_locations.push_back(
                std::make_pair(patch_idx, breeder_idx)
            );
        }
    }

    // Randomise the order so that the first requested locations
    // form a random sample without replacement.
    std::shuffle(
        individual_locations.begin(),
        individual_locations.end(),
        rng_r
    );

    unsigned n_to_sample = param.assay_n_individuals;

    if (n_to_sample > individual_locations.size())
    {
        std::cerr
            << "Requested "
            << n_to_sample
            << " assay individuals, but population contains only "
            << individual_locations.size()
            << ". Sampling the whole population instead.\n";

        n_to_sample =
            static_cast<unsigned>(individual_locations.size());
    }

    // Separate files for individual metadata and hormone trajectories.
    std::ofstream metadata_file{
        param.file_name + "_assay_metadata"
    };

    std::ofstream trajectory_file{
        param.file_name + "_assay_trajectories"
    };

    if (!metadata_file)
    {
        std::cerr
            << "Could not open assay metadata file: "
            << param.file_name + "_assay_metadata"
            << "\n";

        return;
    }

    if (!trajectory_file)
    {
        std::cerr
            << "Could not open assay trajectory file: "
            << param.file_name + "_assay_trajectories"
            << "\n";

        return;
    }

    metadata_file
        << "sample_id;"
        << "patch_index;"
        << "breeder_index;"
        << "baseline_influx0;"
        << "baseline_influx1;"
        << "h1_S0;"
        << "h1_S1;"
        << "removal0;"
        << "removal1;"
        << "baseline_phenotype;"
        << "h1_S_phenotype;"
        << "removal_phenotype;"
        << "equilibrium_hormone"
        << '\n';

    trajectory_file
        << "sample_id;"
        << "patch_index;"
        << "breeder_index;"
        << "relative_time;"
        << "attacked;"
        << "stress_hormone"
        << '\n';

    for (unsigned sample_idx = 0;
         sample_idx < n_to_sample;
         ++sample_idx)
    {
        unsigned patch_idx =
            individual_locations[sample_idx].first;

        unsigned breeder_idx =
            individual_locations[sample_idx].second;

        const Individual &sampled_individual =
            metapopulation[patch_idx].breeders[breeder_idx];

        // TABORSKY VALIDATION:
        // Express diploid traits as the mean of the two alleles,
        // matching the phenotype convention used in the main model.
        double baseline_phenotype =
            0.5 * (sampled_individual.baseline_influx[0] +
                   sampled_individual.baseline_influx[1]);
        
        double h1_S_phenotype =
            0.5 * (sampled_individual.h1_S[0] +
                   sampled_individual.h1_S[1]);
        
        double removal_phenotype =
            0.5 * (sampled_individual.removal[0] +
                   sampled_individual.removal[1]);

        double hormone;

        if (removal_phenotype > 0.0)
        {
            hormone = baseline_phenotype / removal_phenotype;
        }
        
        else
        {
            hormone = param.hmax;
        }

        // Match the hormone clipping used in the main model.
        if (hormone < 0.0)
        {
            hormone = 0.0;
        }

        if (hormone > param.hmax)
        {
            hormone = param.hmax;
        }

        unsigned sample_id = sample_idx + 1;

        metadata_file
            << sample_id << ";"
            << patch_idx << ";"
            << breeder_idx << ";"
            << sampled_individual.baseline_influx[0] << ";"
            << sampled_individual.baseline_influx[1] << ";"
            << sampled_individual.h1_S[0] << ";"
            << sampled_individual.h1_S[1] << ";"
            << sampled_individual.removal[0] << ";"
            << sampled_individual.removal[1] << ";"
            << baseline_phenotype << ";"
            << h1_S_phenotype << ";"
            << removal_phenotype << ";"
            << hormone
            << '\n';

        int first_time =
            -static_cast<int>(param.assay_pre_time);

        int final_time =
            static_cast<int>(param.assay_post_time);

        // State variables for the standardised stress response.
        // The assay begins outside an active response, with no previous
        // stress-induced peak.
        double hx = 0.0;
        unsigned time_since_last_stressor =
            param.tmax_stress_influx;

        // Record the equilibrium baseline throughout the pre-attack period.
        for (int relative_time = first_time;
             relative_time <= final_time;
             ++relative_time)
        {
            bool attacked = relative_time == 0;

            // During the pre-attack period, hormone remains at its
            // equilibrium baseline. At t = 0 a standardised attack
            // starts the Taborsky-style prolonged response.
            if (relative_time >= 0)
            {
                if (attacked)
                {
                    time_since_last_stressor = 0;
                }
                else
                {
                    ++time_since_last_stressor;
                }
            
                // Baseline hormone dynamics.
                hormone =
                    (1.0 - removal_phenotype) * hormone +
                    baseline_phenotype;
            
                // TABORSKY VALIDATION:
                // Stress-induced influx can continue for tmax_stress_influx
                // timesteps following the attack and is regulated by h1_S.
                if (time_since_last_stressor <
                    param.tmax_stress_influx)
                {
                    hormone +=
                        param.stress_influx_max *
                        stress_feedback(hx, h1_S_phenotype);
                }
                }

            if (hormone < 0.0)
            {
                hormone = 0.0;
            }

            if (hormone > param.hmax)
            {
                hormone = param.hmax;
            }
            
            // Track the highest hormone level reached during the
            // standardised response for negative feedback.
            if (relative_time >= 0 && hx < hormone)
            {
                hx = hormone;
            }

            trajectory_file
                << sample_id << ";"
                << patch_idx << ";"
                << breeder_idx << ";"
                << relative_time << ";"
                << attacked << ";"
                << hormone
                << '\n';
        }
    }
}
 
void StressSocial::write_data_headers()

{
	data_file << "time;meanv;varv;"
            << "mean_baseline_influx;var_baseline_influx;"
            << "mean_stress_influx;var_stress_influx;"
            << "mean_h1_S;var_h1_S;"
            << "mean_hstart;var_hstart;"
            << "mean_vigilance_influx;var_vigilance_influx;"
            << "mean_removal;var_removal;"
            << "mean_damage;var_damage;"
            << "mean_stress_hormone;var_stress_hormone;"
            << "total_global_fecundity;"
            << "predator_presence_fraction;" // addition of predator presence fraction in output
            << "mean_group_V;" // mean group vigilance
            << "n_attacked;n_death_damage;mean_damage_at_damage_death;"
            << "n_death_predator;attack_death_fraction;ntotalalive"  
            << std::endl;
            
            // Note: mean_vigilance column is the expressed vigilance phenotype
            // effective_vigilance() = 0.5*(v0+v1) clamped to [0,1], not raw sum of alleles
            
}	

// means and the variances of the various traits
// both the genetic traits and also the non-genetic
// traits. 
void StressSocial::write_data() 
{
    // allocate variables that contain the means
    // allocate variables that contain the sum of squares (for the variances)
    // then calculate the variance as var(x) = sum_of_squares/n - mean(x) * mean(x)

    // Allocate variables that contain means, ss and variance
    int total_individuals{0}; // total number of indiv counter
    double meanv {0.0}; // mean of vigilance
    double ssv {0.0}; // sum of squares of vigilance
    double varv {0.0}; // variance in vigilance
    double mean_baseline_influx {0.0}; // mean baseline influx
    double ss_baseline_influx {0.0}; // sum of squares baseline influx
    double var_baseline_influx {0.0}; // variance in baseline influx
    double mean_stress_influx {0.0}; // mean stress influx
    double ss_stress_influx {0.0}; // sum of squares stress influx
    double var_stress_influx {0.0}; // variance in stress influx
    // TABORSKY VALIDATION:
    // Summary statistics for the evolved h1_S feedback trait.
    double mean_h1_S {0.0};
    double ss_h1_S {0.0};
    double var_h1_S {0.0};
    double mean_hstart{0.0};
    double ss_hstart{0.0};
    double var_hstart{0.0};
    double mean_vigilance_influx {0.0}; // mean vigilance influx
    double ss_vigilance_influx {0.0}; // sum of squares vigilance influx
    double var_vigilance_influx {0.0}; // variance in vigilance influx
    double mean_removal {0.0}; // mean stress hormone removal
    double ss_removal {0.0}; // sum of squares stress hormone removal
    double var_removal {0.0}; // variance in stress hormone removal
    double mean_damage {0.0}; // mean damage
    double ss_damage {0.0}; // sum of squares damage
    double var_damage {0.0}; // variance in damage
    double mean_stress_hormone {0.0}; // mean stress hormone
    double ss_stress_hormone {0.0}; // sum of squares stress hormone
    double var_stress_hormone {0.0}; // variance in stress hormone
    double predator_presence_fraction {0.0}; // track predator presence across patches
    double mean_group_V {0.0}; // mean group vigilance

    for (auto &patch : metapopulation) {
        for (auto &breeder : patch.breeders) {

    // record diploid trait values as average of the two alleles
    // use expressed vigilance phenotype (0.5*(v0+v1, clamped)
            double vigilance = effective_vigilance(breeder, param.vigilance); 
            double baseline_influx = 0.5 * (breeder.baseline_influx[0] + breeder.baseline_influx[1]);
            double stress_influx = 0.5 * (breeder.stress_influx[0] + breeder.stress_influx[1]);
            double vigilance_influx = 0.5 * (breeder.vigilance_influx[0] + breeder.vigilance_influx[1]);
            double h1_S = 0.5 * (breeder.h1_S[0] + breeder.h1_S[1]);
            
            // TABORSKY VALIDATION:
            // Expressed starting hormone phenotype.
            double hstart = 0.5 * (breeder.hstart[0] + breeder.hstart[1]);
            
            double removal = 0.5 * (breeder.removal[0] + breeder.removal[1]);
            double damage = breeder.damage;
            double stress_hormone = breeder.stress_hormone;

            meanv += vigilance;
            ssv += vigilance * vigilance;

            mean_baseline_influx += baseline_influx;
            ss_baseline_influx += baseline_influx * baseline_influx;

            mean_stress_influx += stress_influx;
            ss_stress_influx += stress_influx * stress_influx;

            mean_h1_S += h1_S;
            ss_h1_S += h1_S * h1_S;
            
            mean_hstart += hstart;
            ss_hstart += hstart * hstart;
            
            mean_vigilance_influx += vigilance_influx;
            ss_vigilance_influx += vigilance_influx * vigilance_influx; // Is initialised correctly in individual.cpp?

            mean_removal += removal;
            ss_removal += removal * removal;

            mean_damage += damage;
            ss_damage += damage * damage;

            mean_stress_hormone += stress_hormone;
            ss_stress_hormone += stress_hormone * stress_hormone; 

            ++total_individuals;
        }
    }
    
    // fraction of patches with predator present this timestep - needed in output
    int predator_patches = 0;
    for (const auto &patch : metapopulation) {
        if (patch.predator_patch) ++predator_patches;
    }
    predator_presence_fraction =
        static_cast<double>(predator_patches) / param.npatches;
        
    // mean group vigilance across patches
    for (const auto &patch : metapopulation) {
        mean_group_V += patch.V;
    }

    mean_group_V /= param.npatches;

    // Calculate mean
        if (total_individuals > 0) {
        meanv /= total_individuals;
        mean_baseline_influx /= total_individuals;
        mean_stress_influx /= total_individuals;
        mean_h1_S /= total_individuals;
        mean_hstart /= total_individuals;
        mean_vigilance_influx /= total_individuals;
        mean_removal /= total_individuals;
        mean_damage /= total_individuals;
        mean_stress_hormone /= total_individuals;
    }

     // Calculate variance if total_individuals is not zero
    varv = (total_individuals > 0) ? (ssv / total_individuals - meanv * meanv) : 0.0;
    var_baseline_influx = (total_individuals > 0) ? (ss_baseline_influx / total_individuals - mean_baseline_influx * mean_baseline_influx): 0.0;
    var_stress_influx = (total_individuals > 0) ? (ss_stress_influx / total_individuals - mean_stress_influx * mean_stress_influx) : 0.0;
    var_h1_S = (total_individuals > 0) ? (ss_h1_S / total_individuals - mean_h1_S * mean_h1_S) : 0.0;
    var_hstart = (total_individuals > 0) ? (ss_hstart / total_individuals - mean_hstart * mean_hstart) : 0.0;
    var_vigilance_influx = (total_individuals > 0) ? (ss_vigilance_influx / total_individuals - mean_vigilance_influx * mean_vigilance_influx) : 0.0;
    var_removal = (total_individuals > 0) ? (ss_removal / total_individuals - mean_removal * mean_removal) : 0.0;
    var_damage = (total_individuals > 0) ? (ss_damage / total_individuals - mean_damage * mean_damage) : 0.0;
    var_stress_hormone = (total_individuals > 0) ? (ss_stress_hormone/ total_individuals - mean_stress_hormone * mean_stress_hormone) : 0.0;

    double attack_death_fraction = 0.0;
    
    if (n_attacked > 0)
    {
        attack_death_fraction = static_cast<double>(n_death_predator) / n_attacked;
    }
    
    // Mean damage among individuals that died from damage this timestep.
    // If no individuals died from damage, leave as 0.0.
    double mean_damage_at_damage_death = 0.0;
    
    if (n_death_damage > 0)
    {
        mean_damage_at_damage_death =
            sum_damage_at_damage_death / n_death_damage;
}

    unsigned int ntotal = param.npatches * param.n;

    assert(ntotal >= n_death_damage + n_death_predator); //error checking as ntotal should always be greater
    
    data_file << time_step << ";"
        << meanv << ";" 
        << varv << ";" 
        << mean_baseline_influx << ";"
        << var_baseline_influx << ";" 
        << mean_stress_influx << ";"
        << var_stress_influx << ";" 
        << mean_h1_S << ";"
        << var_h1_S << ";"
        << mean_hstart << ";"
        << var_hstart << ";"
        << mean_vigilance_influx << ";"
        << var_vigilance_influx << ";" 
        << mean_removal << ";"
        << var_removal << ";" 
        << mean_damage << ";"
        << var_damage << ";" 
        << mean_stress_hormone << ";"
        << var_stress_hormone << ";" 
        << last_total_global_fecundity << ";"
        << predator_presence_fraction << ";"
        << mean_group_V << ";"
        << n_attacked << ";"
        << n_death_damage << ";"
        << mean_damage_at_damage_death << ";"
        << n_death_predator << ";"
        << attack_death_fraction << ";"
        << (ntotal - n_death_damage - n_death_predator)
        << '\n';
}


void StressSocial::predator_visit()
{
    double V; // auxiliary variable reflecting whether at least a single individual is vigilant

    unsigned random_breeder_idx;

    // 1. all patches that are of type P need to have a visit by a predator
    // 2. predator samples x individuals to attack 
    // 3. predator attacks them, so this changes an individuals' state
    // 4. an individual can avoid attack dependent on its strrress response
    for (auto metapop_iter = metapopulation.begin();
            metapop_iter != metapopulation.end();
            ++metapop_iter)
    {
        if (metapop_iter->predator_patch)
        {
            // check whether at least one individual is vigilant
            V = calculate_group_vigilance(*metapop_iter);

            // calculate the probability that nobody is vigilant - attacks happen when no group members are vigilant
            // higher value of V = fewer attacks
            if (uniform(rng_r) < 1.0 - V && 
                    uniform(rng_r) < param.p_attack)
            {
                // then sample which individual will die
                random_breeder_idx = take_random_breeder(rng_r);

                metapop_iter->breeders[random_breeder_idx].is_attacked = true;
                
                // TABORSKY VALIDATION:
                // An attack starts/restarts the prolonged stress-response window.
                metapop_iter->breeders[random_breeder_idx].time_since_last_stressor = 0;

                ++n_attacked;
                
                // debug: checking why all attacked individuals die
                if (time_step == 100 ||
                    time_step == 1000 ||
                    (time_step >= 10000 && time_step % 10000 == 0))
                {                
                std::cout << time_step
                          << " h "
                          << metapop_iter->breeders[random_breeder_idx].stress_hormone
                          << " attack_survival "
                          << attack_survival(
                                metapop_iter->breeders[random_breeder_idx].stress_hormone)
                          << std::endl;
                }
   
                if (
                        uniform(rng_r) < 
                            attack_survival(metapop_iter->breeders[random_breeder_idx].stress_hormone))
                {
                    // we need to implement that individuals can flee the attack 
                    // dependent on their stress hormone level h
                    metapop_iter->breeders[random_breeder_idx].is_alive = true;
                }
                else
                {
                    ++n_death_predator;
                    metapop_iter->breeders[random_breeder_idx].is_alive = false;
                }
            }
            
            
            // then store V in the patch object
            metapop_iter->V = V;
        }
        else
        {
            // patch has no predator this timestep
            // set V explicitly to 0.0 so it’s always initialised
            metapop_iter->V = 0.0;
        }
    }
} // end predator_visit()



// Probability of surviving an attack given hormone level h,
// With hmax = 1 and survival_power = 1 this is equivalent to
// the Taborsky stress model
double StressSocial::attack_survival(double const h)
{
    return(pow(h/param.hmax, param.survival_power));
}

// TABORSKY VALIDATION:
// Hormone-dependent negative feedback on stress-induced influx.
// If previous peak hormone hx exceeds h1_S, stress-induced production stops.
// Otherwise production declines linearly as hx approaches h1_S.
double StressSocial::stress_feedback(
        double const hx,
        double const h1_S)
{
    if (h1_S <= 0.0 || hx > h1_S)
    {
        return 0.0;
    }

    return 1.0 - hx / h1_S;
}

// go over all patches and calculate the total probability
// that none of the individuals are vigilant. 
double StressSocial::calculate_group_vigilance(Patch const &current_patch)
{
    double prob_none_vigilant = 1.0;

    // go over all patches and calculate the total probability
    // that none of the individuals are vigilant
    for (auto breeder_iter = current_patch.breeders.begin();
            breeder_iter != current_patch.breeders.end();
            ++breeder_iter)
    {

    // Use expressed vigilance phenotype (bounded [0,1])
        double v_eff = effective_vigilance(*breeder_iter, param.vigilance);
        prob_none_vigilant = prob_none_vigilant * (1.0 - v_eff);
    }

    // the probability that none of the individuals 
    // are vigilant is now calculated after this loop.
    //
    // from this, we can then get the probability that
    // at least 1 individual is vigilant

    // return 1 - (1-v)^n
    return 1.0 - prob_none_vigilant;
} // end calculate_group_vigilance()


// calculate how damage affects survival
void StressSocial::survive_damage_vigilance()
{
    double d;

    for (auto metapop_iter = metapopulation.begin();
            metapop_iter != metapopulation.end();
            ++metapop_iter)
    {
        // loop over all breeders, evaluate damage
        // kill them if they die
        for (unsigned breeder_idx{0};
                breeder_idx < metapop_iter->breeders.size();
                ++breeder_idx)
        {
            if (metapop_iter->breeders[breeder_idx].is_alive)
            {
                    d = metapop_iter->breeders[breeder_idx].damage;

                    assert(std::isfinite(d));

            // use expressed vigilance phenotype (bounded [0,1])
                    double v = effective_vigilance(metapop_iter->breeders[breeder_idx], param.vigilance);


                    if (uniform(rng_r) < 1.0 - mu(d, v))
                    {
                        // note that mortality due to lack of vigilance
                        // is elsewhere, namely in predator_visit()
                        // individual does not survive
                metapop_iter->breeders[breeder_idx].is_alive = true;
                    }
                    else
                    {
                      if (time_step == 100 ||
                          time_step == 1000 ||
                          (time_step >= 10000 && time_step % 10000 == 0))
                      {
                          std::cout << time_step
                                    << " DAMAGE_DEATH "
                                    << "damage "
                                    << d
                                    << " vigilance "
                                    << v
                                    << " attacked "
                                    << metapop_iter->breeders[breeder_idx].is_attacked
                                    << std::endl;
                      }
                        
                        // Record the individual's damage immediately before it dies from
                        // damage-related mortality. This is accumulated across all such deaths
                        // during the timestep and converted to a mean in write_data().
                        sum_damage_at_damage_death += d;
                                  
                        ++n_death_damage;
                    metapop_iter->breeders[breeder_idx].is_alive = false;
                }
            } // end if metapop_iter
        } // end for unsigned breeder_idx
    } // for end metapop_iter 
} // end survival_damage()


double StressSocial::mu(
        double const damage,
        double const vigilance)
{
    // TABORSKY VALIDATION:
    // Damage no longer contributes directly to mortality.
    // In the Taborsky stress model, elevated hormone/damage reduces
    // reproductive success rather than increasing mortality.
    //
    // Keep the damage argument for now because it is part of the existing
    // function interface and may be restored when this validation branch
    // is compared with the full social model.
    (void)damage;

    double mortality_prob{
        param.m0 + param.mv * vigilance
    };

    return(mortality_prob);
}


void StressSocial::reproduce()
{
    // list of all the group level fecundities
    std::vector <double> group_level_fecundities;

    // total fecundity across all groups
    double total_global_fecundity = 0.0;

    // list of all the fecundities across all breeders of a single patch
    std::vector <double> individual_level_fecundities;

    // auxiliary variable to calculate group level fecundity
    double group_level_fecundity, individual_fecundity;

    // calculate a mean fecundity distribution
    for (auto metapop_iter = metapopulation.begin();
            metapop_iter != metapopulation.end();
            ++metapop_iter)
    {
        // reset the group level total fecundity
        // as we start with a new patch
        group_level_fecundity = 0.0;

        // reset the vector with individual level fecundities
        // as we start with a new patch
        individual_level_fecundities.clear();

        // calculate fecundity for each group
        // dependent on individual vigilance values
        for (auto breeder_iter = metapop_iter->breeders.begin();
                breeder_iter != metapop_iter->breeders.end();
                ++breeder_iter)
        {
        
            // TABORSKY VALIDATION:
            // Dead individuals must have zero reproductive success.
            // In the Taborsky stress model, dead individuals remain in the
            // population until replacement but have fecundity exactly zero.
            if (!breeder_iter->is_alive)
            {
                individual_fecundity = 0.0;
            }
            else
            {
                // Expressed vigilance phenotype, bounded to [0,1].
                double v_eff =
                    effective_vigilance(*breeder_iter, param.vigilance);
            
                // Existing fecundity cost of vigilance.
                double vigilance_fecundity =
                    1.0 - std::pow(
                        v_eff,
                        param.fecundity_power
                    );
            
                // Taborsky damage-dependent fecundity:
                // F = 1 - (damage / dmax)^ad
                double damage_fecundity =
                    1.0 - std::pow(
                        breeder_iter->damage / param.dmax,
                        param.damage_fecundity_power
                    );
            
                damage_fecundity =
                    std::clamp(
                        damage_fecundity,
                        0.0,
                        1.0
                    );
            
                individual_fecundity =
                    vigilance_fecundity *
                    damage_fecundity;
            }


            individual_level_fecundities.push_back(individual_fecundity);

            group_level_fecundity += individual_fecundity;
        }

        // param object to update this patch's fecundity distribution
        std::discrete_distribution<unsigned>::param_type 
            fecundity_distribution_param(
                    individual_level_fecundities.begin()
                    ,individual_level_fecundities.end());

        // update the patch's discrete distribution of fecundities
        // with this param_type object we just made
        metapop_iter->within_patch_fecundity_distribution.param(
                fecundity_distribution_param);

        // add the total fecundity value to the list of group-level fecundities
        group_level_fecundities.push_back(group_level_fecundity);
        

        // accumulate global fecundity across all patches so that migration works
        total_global_fecundity += group_level_fecundity;

    }


    // make a probability distribution of the patch level fecundities
    std::discrete_distribution<unsigned> group_level_fecundity_distribution(
            group_level_fecundities.begin(), 
            group_level_fecundities.end());

    // variable holding the patch we will sample new offspring from
    unsigned patch_producing_new_offspring_idx, mum_idx, dad_idx;
    double probability_sample_immigrant, migrant_contribution, local_contribution, total_local_fecundity;

    // TABORSKY VALIDATION:
    // Equilibrium probability that a newly born individual experiences
    // the predator-present environment.
    //
    // Taborsky assigns every newborn a new environmental state drawn
    // independently from this equilibrium probability. In social_stress,
    // environment is a patch property rather than an individual property,
    // so we reproduce this behaviour only when n = 1, where one patch
    // corresponds to one individual.
    double equilibrium_predator_probability =
        param.s[NP] /
        (param.s[NP] + param.s[P]);

    // tasks ahead:
    // 1. go over all breeders
    // 2. are they dead?
    // 3. sample new offspring from distribution
    // 4. replace dead breeder with new offspring.
    // done
    
    // calculate a mean fecundity distribution
    for (unsigned patch_idx = 0;
            patch_idx < param.npatches;
            ++patch_idx)
    {
        total_local_fecundity = group_level_fecundities[patch_idx];

        // total_global_fecundity is now accumulated above
        // so migrant_contribution can be > 0 and immigration can actually occur 
        migrant_contribution = param.p_mig * total_global_fecundity / param.npatches;

        local_contribution = (1.0 - param.p_mig) * total_local_fecundity;

        // probability that we draw parents from another patch
        probability_sample_immigrant = migrant_contribution / (local_contribution + migrant_contribution);

        // calculate fecundity for each group
        // dependent on individual vigilance values
        for (auto breeder_iter = metapopulation[patch_idx].breeders.begin();
                breeder_iter != metapopulation[patch_idx].breeders.end();
                ++breeder_iter)
        {
            if (!breeder_iter->is_alive) // individual dead, hence needs replacing
            {
                // get offspring from remote patch
                if (uniform(rng_r) < probability_sample_immigrant)
                {
                    // sample remote patch
                    patch_producing_new_offspring_idx = group_level_fecundity_distribution(rng_r);

                    assert(patch_producing_new_offspring_idx < param.npatches);

                } else // get offspring from local patch
                {
                    patch_producing_new_offspring_idx = patch_idx;
                    
                    assert(patch_producing_new_offspring_idx < param.npatches);
                }

                // we know the patch, now which parent
                // first pick mum
                mum_idx = metapopulation[
                    patch_producing_new_offspring_idx].within_patch_fecundity_distribution(rng_r);

                assert(mum_idx < param.n);

                // pick dad
                dad_idx = metapopulation[
                    patch_producing_new_offspring_idx].within_patch_fecundity_distribution(rng_r);

                assert(dad_idx < param.n);

                // call birth constructor
                Individual Kid(
                        metapopulation[patch_producing_new_offspring_idx].breeders[mum_idx],
                        metapopulation[patch_producing_new_offspring_idx].breeders[dad_idx],
                        param,
                        rng_r);

                assert(Kid.v[0] >= 0);
                assert(Kid.v[0] <= 1.0);

                // fill the vacancy with new offspring
                *breeder_iter = Kid;
                
                // TABORSKY VALIDATION:
                // In Taborsky, each newborn starts in a newly sampled environmental
                // state rather than inheriting the environmental history of the
                // individual it replaces.
                //
                // social_stress stores predator state at patch level, so this is
                // equivalent only for the validation setup n = 1. For social-model
                // runs with n > 1, retain the existing shared patch environment.
                if (param.n == 1)
                {
                    metapopulation[patch_idx].predator_patch =
                        uniform(rng_r) <
                        equilibrium_predator_probability;
                }
                
                assert(breeder_iter->v[0] >= 0);
                assert(breeder_iter->v[0] <= 1.0);

            } // end if breeder_iter is_alive
        } // end for breeder_iter
    } // end for patch_idx
    
      // store total fecundity for reporting in write_data()
      last_total_global_fecundity = total_global_fecundity;
    
} // end StressSocial::reproduce()



// write parameters to file
// added in group and population size
void StressSocial::write_parameters() 
{
    data_file << std::endl
        << std::endl
        << "seed;" << seed << ";" << std::endl
        << "time_step;" << time_step << ";" << std::endl
        << "dispersal;" << param.p_mig << ";" << std::endl
        << "npatches;" << param.npatches << ";" << std::endl
        << "n;" << param.n << ";" << std::endl
        << "ntotal;" << param.npatches * param.n << ";" << std::endl
        << "s_np;" << param.s[NP] << ";" << std::endl
        << "s_p;" << param.s[P] << ";" << std::endl
        << "p_attack;" << param.p_attack << ";" << std::endl
        << "fecundity_power;" << param.fecundity_power << ";" << std::endl
        << "damage_fecundity_power;" << param.damage_fecundity_power << ";" << std::endl
        << "hmin;" << param.hmin << ";" << std::endl
        << "hmax;" << param.hmax << ";" << std::endl
        << "dmax;" << param.dmax << ";" << std::endl
        << "survival_power;" << param.survival_power << ";" << std::endl
        << "vigilance;" << param.vigilance << ";" << std::endl // vigilance/on off written to file (0 off, 1 on)
        << "run_end_assay;" << param.run_end_assay << ";" << std::endl
        << "assay_n_individuals;" << param.assay_n_individuals << ";" << std::endl
        << "assay_pre_time;" << param.assay_pre_time << ";" << std::endl
        << "assay_post_time;" << param.assay_post_time << ";" << std::endl
        << "init_v;" << param.init_v << ";" << std::endl
        << "init_stress_hormone_level;" << param.init_stress_hormone_level << ";" << std::endl // A legacy in this model
        << "init_baseline_influx;" << param.init_baseline_influx << ";" << std::endl
        << "init_removal;" << param.init_removal << ";" << std::endl
        << "min_removal;" << param.min_removal << ";" << std::endl
        << "init_hstart;" << param.init_hstart << ";" << std::endl
        << "mu_baseline;" << param.mu_baseline << ";" << std::endl
        << "mu_stress_influx;" << param.mu_stress_influx << ";" << std::endl
        << "mu_vigilance_influx;" << param.mu_vigilance_influx << ";" << std::endl
        << "mu_removal;" << param.mu_removal << ";" << std::endl
        << "mu_v;" << param.mu_v << ";" << std::endl
        << "mu_hstart;" << param.mu_hstart << ";" << std::endl
        << "sdmu;" << param.sdmu << ";" << std::endl
        << "stress_influx_max;" << param.stress_influx_max << ";" << std::endl
        << "tmax_stress_influx;" << param.tmax_stress_influx << ";" << std::endl
        << "init_h1_S;" << param.init_h1_S << ";" << std::endl
        << "mu_h1_S;" << param.mu_h1_S << ";" << std::endl
        << "file_name;" << param.file_name << ";" << std::endl 
        << std::endl
        << std::endl;
}

// TABORSKY VALIDATION:
// First physiological phase of each timestep.
//
// This matches the beginning of Taborsky's survive() function:
// 1. advance time since the previous stressor;
// 2. apply normal hormone removal and baseline influx.
//
// This occurs BEFORE the predator attack, so attack survival uses
// the individual's hormone level after the current baseline update.
void StressSocial::update_baseline_hormone()
{
    for (auto metapop_iter = metapopulation.begin();
         metapop_iter != metapopulation.end();
         ++metapop_iter)
    {
        for (auto breeder_iter = metapop_iter->breeders.begin();
             breeder_iter != metapop_iter->breeders.end();
             ++breeder_iter)
        {
            // All vacancies were replaced at the end of the previous
            // timestep, so individuals should normally be alive here.
            if (!breeder_iter->is_alive)
            {
                continue;
            }

            // Taborsky increments this counter at the beginning
            // of each individual's timestep, before any new attack.
            ++breeder_iter->time_since_last_stressor;

            double removal =
                0.5 * (
                    breeder_iter->removal[0] +
                    breeder_iter->removal[1]
                );

            double baseline_influx =
                0.5 * (
                    breeder_iter->baseline_influx[0] +
                    breeder_iter->baseline_influx[1]
                );

            double vigilance_influx =
                0.5 * (
                    breeder_iter->vigilance_influx[0] +
                    breeder_iter->vigilance_influx[1]
                );

            // Normal hormone dynamics.
            breeder_iter->stress_hormone =
                (1.0 - removal) *
                breeder_iter->stress_hormone +
                baseline_influx +
                vigilance_influx * metapop_iter->V;

            // Keep hormone within its permitted range.
            breeder_iter->stress_hormone =
                std::clamp(
                    breeder_iter->stress_hormone,
                    0.0,
                    param.hmax
                );
        }
    }
}

// TABORSKY VALIDATION:
// Second physiological phase of each timestep.
//
// predator_visit() has already happened at this point.
// Any attacked individual therefore has
// time_since_last_stressor = 0.
//
// Stress-induced hormone production then continues while the
// individual remains inside the post-stressor response window.
void StressSocial::update_stress_response()
{
    for (auto metapop_iter = metapopulation.begin();
         metapop_iter != metapopulation.end();
         ++metapop_iter)
    {
        for (auto breeder_iter = metapop_iter->breeders.begin();
             breeder_iter != metapop_iter->breeders.end();
             ++breeder_iter)
        {
            // Individuals killed by the attack no longer contribute
            // to later survival, damage or reproduction. There is no
            // need to update their physiological state further.
            if (!breeder_iter->is_alive)
            {
                breeder_iter->is_attacked = false;
                continue;
            }

            if (breeder_iter->time_since_last_stressor <
                param.tmax_stress_influx)
            {
                double h1_S =
                    0.5 * (
                        breeder_iter->h1_S[0] +
                        breeder_iter->h1_S[1]
                    );

                breeder_iter->stress_hormone +=
                    param.stress_influx_max *
                    stress_feedback(
                        breeder_iter->hx,
                        h1_S
                    );

                breeder_iter->stress_hormone =
                    std::clamp(
                        breeder_iter->stress_hormone,
                        0.0,
                        param.hmax
                    );

                // Match Taborsky: hx is updated only during
                // the active post-stressor response.
                if (breeder_iter->hx <
                    breeder_iter->stress_hormone)
                {
                    breeder_iter->hx =
                        breeder_iter->stress_hormone;
                }
            }

            // Reset attack flag ready for the next timestep.
            // The stress-response timer remains active independently.
            breeder_iter->is_attacked = false;
        }
    }
}

// TABORSKY VALIDATION:
// Final physiological phase before reproduction.
//
// This occurs AFTER predator and background mortality, matching
// Taborsky's lifecycle. Only surviving individuals accumulate
// current-timestep hormone-dependent damage.
//
// d(t+1) = (1-g)d(t) + k*h(t+1)
void StressSocial::update_damage()
{
    for (auto metapop_iter = metapopulation.begin();
         metapop_iter != metapopulation.end();
         ++metapop_iter)
    {
        for (auto breeder_iter = metapop_iter->breeders.begin();
             breeder_iter != metapop_iter->breeders.end();
             ++breeder_iter)
        {
            // Taborsky updates damage only for individuals
            // that remain alive after mortality.
            if (!breeder_iter->is_alive)
            {
                continue;
            }

            double damage_tplus1 =
                (1.0 - param.g) *
                breeder_iter->damage +
                param.k *
                breeder_iter->stress_hormone;

            breeder_iter->damage =
                std::clamp(
                    damage_tplus1,
                    0.0,
                    param.dmax
                );
        }
    }
}
