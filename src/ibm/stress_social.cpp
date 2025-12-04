// heart of the stress social code

#include <cassert>
#include <iostream>
#include <cmath>
#include <numeric>
#include <fstream>

#include "stress_social.hpp"
#include "patch.hpp"
#include "individual.hpp"

inline double effective_vigilance(Individual const & ind)
{
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
    // EG added writing of distribution to new file at time_step == 0 
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

          // effectively, we want the predator to visit some patches
          // and to attack individuals there.
            
            switch_predator_status();

            predator_visit();
            survive_damage_vigilance();
            reproduce();
            update_stress_hormone();
        
           // error checking: ntotal should always be >= each death count - simplified from previous version
            assert(param.n * param.npatches >= n_death_damage); // total pop >= deaths from damage
            assert(param.n * param.npatches >= n_death_predator); // total pop >= deaths from predation
            assert(param.n * param.npatches >= n_death_damage + n_death_predator); // total pop >= total deaths
        
            if (time_step % param.data_output_interval == 0)
            {
                write_data();
                // write_distribution(); // debug only
            }
    } // end for time_step
    
    write_parameters();
} // end StressSocial constructor

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

// EG FIX: update whether each patch has a predator, once per timestep
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
// EG added: creation of separate debug file
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
 
void StressSocial::write_data_headers()

{
	data_file << "time;seed;meanv;varv;"
            << "mean_baseline_influx;var_baseline_influx;"
            << "mean_stress_influx;var_stress_influx;"
            << "mean_vigilance_influx;var_vigilance_influx;"
            << "mean_removal;var_removal;"
            << "mean_damage;var_damage;"
            << "mean_stress_hormone;var_stress_hormone;"
            << "total_global_fecundity;"
            << "predator_presence_fraction;" // addition of predator presence fraction in output
            << "n_attacked;n_death_damage;n_death_predator;ntotalalive" 
            << std::endl;
            
            // EG NOTE: mean_vigilance column is the expressed vigilance phenotype
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
    double predator_presence_fraction = 0.0; // track predator presence across patches

    for (auto &patch : metapopulation) {
        for (auto &breeder : patch.breeders) {

    // EG FIX: record diploid trait values as average of the two alleles
    // EG FIX: Use expressed vigilance phenotype (0.5*(v0+v1, clamped)
            double vigilance = effective_vigilance(breeder); 
            double baseline_influx = 0.5 * (breeder.baseline_influx[0] + breeder.baseline_influx[1]);
            double stress_influx = 0.5 * (breeder.stress_influx[0] + breeder.stress_influx[1]);
            double vigilance_influx = 0.5 * (breeder.vigilance_influx[0] + breeder.vigilance_influx[1]);
            double removal = 0.5 * (breeder.removal[0] + breeder.removal[1]);
            double damage = breeder.damage;
            double stress_hormone = breeder.stress_hormone;

            meanv += vigilance;
            ssv += vigilance * vigilance;

            mean_baseline_influx += baseline_influx;
            ss_baseline_influx += baseline_influx * baseline_influx;

            mean_stress_influx += stress_influx;
            ss_stress_influx += stress_influx * stress_influx;

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
    
    // EG FIX: Fraction of patches with predator present this timestep - needed in output
    int predator_patches = 0;
    for (const auto &patch : metapopulation) {
        if (patch.predator_patch) ++predator_patches;
    }
    predator_presence_fraction =
        static_cast<double>(predator_patches) / param.npatches;

    // Calculate mean
    
        if (total_individuals > 0) {
        meanv /= total_individuals;
        mean_baseline_influx /= total_individuals;
        mean_stress_influx /= total_individuals;
        mean_vigilance_influx /= total_individuals;
        mean_removal /= total_individuals;
        mean_damage /= total_individuals;
        mean_stress_hormone /= total_individuals;
    }

     // Calculate variance if total_individuals is not zero
    varv = (total_individuals > 0) ? (ssv / total_individuals - meanv * meanv) : 0.0;
    var_baseline_influx = (total_individuals > 0) ? (ss_baseline_influx / total_individuals - mean_baseline_influx * mean_baseline_influx): 0.0;
    var_stress_influx = (total_individuals > 0) ? (ss_stress_influx / total_individuals - mean_stress_influx * mean_stress_influx) : 0.0;
    var_vigilance_influx = (total_individuals > 0) ? (ss_vigilance_influx / total_individuals - mean_vigilance_influx * mean_vigilance_influx) : 0.0;
    var_removal = (total_individuals > 0) ? (ss_removal / total_individuals - mean_removal * mean_removal) : 0.0;
    var_damage = (total_individuals > 0) ? (ss_damage / total_individuals - mean_damage * mean_damage) : 0.0;
    var_stress_hormone = (total_individuals > 0) ? (ss_stress_hormone/ total_individuals - mean_stress_hormone * mean_stress_hormone) : 0.0;

    unsigned int ntotal = param.npatches * param.n;

    assert(ntotal >= n_death_damage + n_death_predator); //error checking as ntotal should always be greater
    
    data_file << time_step << ";"
        << seed << ";"
        << meanv << ";" 
        << varv << ";" 
        << mean_baseline_influx << ";"
        << var_baseline_influx << ";" 
        << mean_stress_influx << ";"
        << var_stress_influx << ";" 
        << mean_vigilance_influx << ";"
        << var_vigilance_influx << ";" 
        << mean_removal << ";"
        << var_removal << ";" 
        << mean_damage << ";"
        << var_damage << ";" 
        << mean_stress_hormone << ";"
        << var_stress_hormone << ";" 
        << last_total_global_fecundity << ";"
        << predator_presence_fraction << ";" // addition
        << n_attacked << ";"
        << n_death_damage << ";"
        << n_death_predator << ";"
        << (ntotal - n_death_damage - n_death_predator)
        << '\n'; // EG fix issue with semicolon stuck after death info 

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

            // calculate the probability that nobody is vigilant
            if (uniform(rng_r) < 1.0 - V && 
                    uniform(rng_r) < param.p_attack)
            {
                // then sample which individual will die
                random_breeder_idx = take_random_breeder(rng_r);

                metapop_iter->breeders[random_breeder_idx].is_attacked = true;

                ++n_attacked;

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
            // EG FIX: patch has no predator this timestep
            // set V explicitly to 0.0 so it’s always initialised
            metapop_iter->V = 0.0;
        }
    }
} // end predator_visit()



// probability of surviving an attack given hormone level h
double StressSocial::attack_survival(double const h)
{
    return(pow(h/param.hmax, param.survival_power));
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

    // EG FIX: Use expressed vigilance phenotype (bounded [0,1])
        double v_eff = effective_vigilance(*breeder_iter);
        prob_none_vigilant = prob_none_vigilant * (1.0 - v_eff);
    }

    // ok, now the probability that none of the individuals 
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

            // EG FIX: use expressed vigilance phenotype (bounded [0,1])
                    double v = effective_vigilance(metapop_iter->breeders[breeder_idx]);

                    if (uniform(rng_r) < 1.0 - mu(d, v))
                    {
                        // note that mortality due to lack of vigilance
                        // is elsewhere, namely in predator_visit()
                        // individual does not survive
                metapop_iter->breeders[breeder_idx].is_alive = true;
                    }
                    else
                    {
                        ++n_death_damage;
                    metapop_iter->breeders[breeder_idx].is_alive = false;
                }
            } // end if metapop_iter
        } // end for unsigned breeder_idx
    } // for end metapop_iter 
} // end survival_damage()


// mortality due to baseline mortality, damage and investment in vigilance
double StressSocial::mu(
        double const damage,
        double const vigilance)
{
    double mortality_prob{param.m0 + param.md * damage / param.dmax + param.mv * vigilance};

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
            // calculate 1 - v^x

            double v_eff = effective_vigilance(*breeder_iter);
    // EG FIX: fecundity cost uses expressed vigilance phenotype (bounded [0,1])
            individual_fecundity = 1.0 - std::pow(v_eff, param.fecundity_power);


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
        

        // FIX EG: accumulate global fecundity across all patches so that migration works
        total_global_fecundity += group_level_fecundity;

    }


    // make a probability distribution of the patch level fecundities
    std::discrete_distribution<unsigned> group_level_fecundity_distribution(
            group_level_fecundities.begin(), 
            group_level_fecundities.end());

    // variable holding the patch we will sample new offspring from
    unsigned patch_producing_new_offspring_idx, mum_idx, dad_idx;
    double probability_sample_immigrant, migrant_contribution, local_contribution, total_local_fecundity;


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

        // FIX EG: total_global_fecundity is now accumulated above
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
                // first pick mom
                mum_idx = metapopulation[
                    patch_producing_new_offspring_idx].within_patch_fecundity_distribution(rng_r);

                assert(mum_idx < param.n);

                // pick dad
                dad_idx = metapopulation[
                    patch_producing_new_offspring_idx].within_patch_fecundity_distribution(rng_r);

                assert(mum_idx < param.n);

                // call birth constructor
                Individual Kid(
                        metapopulation[patch_producing_new_offspring_idx].breeders[mum_idx],
                        metapopulation[patch_producing_new_offspring_idx].breeders[dad_idx],
                        param,
                        rng_r);

                assert(Kid.v[0] >= 0);
                assert(Kid.v[0] <= 1.0);

                // fill the vacancy with new kid
                *breeder_iter = Kid;
                
                assert(breeder_iter->v[0] >= 0);
                assert(breeder_iter->v[0] <= 1.0);

            } // end if breeder_iter is_alive
        } // end for breeder_iter
    } // end for patch_idx
    
      // EG FIX: store total fecundity for reporting in write_data()
      last_total_global_fecundity = total_global_fecundity;
    
} // end StressSocial::reproduce()



// write parameters to file
void StressSocial::write_parameters() 
{
    data_file << std::endl
        << std::endl
        << "seed;" << seed << ";" << std::endl
        << "time_step;" << time_step << ";" << std::endl
        << "dispersal;" << param.p_mig << ";" << std::endl
        << "npatches;" << param.npatches << ";" << std::endl
        << "s_np;" << param.s[NP] << ";" << std::endl
        << "s_p;" << param.s[P] << ";" << std::endl
        << "p_attack;" << param.p_attack << ";" << std::endl
        << "fecundity_power;" << param.fecundity_power << ";" << std::endl
        << "hmin;" << param.hmin << ";" << std::endl
        << "hmax;" << param.hmax << ";" << std::endl
        << "dmax;" << param.dmax << ";" << std::endl
        << "survival_power;" << param.survival_power << ";" << std::endl
        << "init_v;" << param.init_v << ";" << std::endl
        << "init_stress_hormone_level;" << param.init_stress_hormone_level << ";" << std::endl
        << "init_removal;" << param.init_removal << ";" << std::endl // EG add removal to output
        << "mu_baseline;" << param.mu_baseline << ";" << std::endl
        << "mu_stress_influx;" << param.mu_stress_influx << ";" << std::endl
        << "mu_vigilance_influx;" << param.mu_vigilance_influx << ";" << std::endl
        << "mu_removal;" << param.mu_removal << ";" << std::endl
        << "mu_v;" << param.mu_v << ";" << std::endl
        << "file_name;" << param.file_name << ";" << std::endl 
        << std::endl
        << std::endl;
}

// update the stress hormone level for each individual
void StressSocial::update_stress_hormone()
{
    double stress_hormone_tplus1, stress_hormone, r, stress_influx, 
           vigilance_influx, baseline_influx, damage,damage_tplus1;

    // calculate a mean fecundity distribution
    for (auto metapop_iter = metapopulation.begin();
            metapop_iter != metapopulation.end();
            ++metapop_iter)
    {
        // calculate fecundity for each group
        // dependent on individual vigilance values
        for (auto breeder_iter = metapop_iter->breeders.begin();
                breeder_iter != metapop_iter->breeders.end();
                ++breeder_iter)
        {
            r = breeder_iter->removal[0] + breeder_iter->removal[1];
            baseline_influx = breeder_iter->baseline_influx[0] 
                + breeder_iter->baseline_influx[1];

            stress_influx = breeder_iter->stress_influx[0] 
                + breeder_iter->stress_influx[1];

            vigilance_influx = breeder_iter->vigilance_influx[0] 
                + breeder_iter->vigilance_influx[1];

            damage = breeder_iter->damage;

            stress_hormone = breeder_iter->stress_hormone;

            stress_hormone_tplus1 = (1.0 - r) * stress_hormone + 
                baseline_influx +
                breeder_iter->is_attacked * stress_influx + 
                vigilance_influx * metapop_iter->V;

// TODO EG: add v=f(stress) here

            // d(t+1) = (1-g)d + k(h - theta_h)^2
            damage_tplus1 = (1.0 - param.g) * damage + 
                param.k * (stress_hormone - param.theta_hormone) * (stress_hormone - param.theta_hormone);

            // undo the is_attacked variable, ready for the next time step
            breeder_iter->is_attacked = false;

            breeder_iter->stress_hormone = stress_hormone_tplus1;
            breeder_iter->damage = damage_tplus1;
        }
    }
} // update_stress_hormone()

