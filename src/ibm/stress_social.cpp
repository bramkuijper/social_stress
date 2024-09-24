// heart of the stress social code

#include "stress_social.hpp"

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
{
    // make some patches P and some NP
    initialize_patches();

    // now run the thing
    for (time_step = 0; time_step < param.max_time; ++time_step)
    {
        // do something!
        // effectively, we want the predator to visit some patches
        // and to attack individuals there.
        predator_visit();
        
        survive_damage_vigilance();

        reproduce();

        update_stress_hormone();

        update_damage();
    }
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

// means and the variances of the various traits
// both the genetic traits and also the non-genetic
// traits. 
void StressSocial::write_data() 
{
    // allocate variables that contain the means
    // allocate variables that contain the sum of squares (for the variances)
    // then calculate the variance as var(x) = sum_of_squares/n - mean(x) * mean(x)

    data_file << time_step << std::endl;
}




void StressSocial::predator_visit()
{
    double V; // auxiliary variable reflecting 
              // whether at least a single individual is vigilant
              //
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
            if (uniform(rng_r) < 1.0 - V)
            {
                // then sample which individual will die
                random_breeder_idx = take_random_breeder(rng_r);

                metapop_iter->breeders[random_breeder_idx].is_attacked <- true;

                // we need to implement that individuals can flee the attack 
                // dependent on their stress hormone level h
                // TODO mplement stress hormone dynamic
                if (uniform(rng_r) < 1.0 - attack_survival(metapop_iter->breeders[random_breeder_idx].stress_hormone))
                {
                    // overwrite breeder at position with final breeder in stack
                    metapop_iter->breeders[random_breeder_idx] = metapop_iter->breeders.back();

                    // delete final element
                    metapop_iter->breeders.pop_back();
                }
            }

            // then store V in the patch object
            metapop_iter->V = V;
        }
    }
} // end predator_visit()


// probability of surviving an attack given hormone level h
double StressSocial::attack_survival(double const h)
{
    return(pow(h/param.hmax, param.survival_power));
}

// go over all patches and calculate the total probability
// that none of the individuals are vigilant
double StressSocial::calculate_group_vigilance(Patch const &current_patch)
{
    double prob_none_vigilant = 1.0;

    // go over all patches and calculate the total probability
    // that none of the individuals are vigilant
    for (auto breeder_iter = current_patch.breeders.begin();
            breeder_iter != current_patch.breeders.end();
            ++breeder_iter)
    {
        prob_none_vigilant = prob_none_vigilant * (1.0 - breeder_iter->v[0] - breeder_iter->v[1]);
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
// TODO: we need to update damage levels each and every timestep
void StressSocial::survive_damage_vigilance()
{
    double d,v;

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
            // TODO check compiler warning if you use += instead of =
            d = metapop_iter->breeders[breeder_idx].damage;
            v = metapop_iter->breeders[breeder_idx].v[0] + metapop_iter->breeders[breeder_idx].v[1];

            // TODO think about mortality due to lack of vigilance
            // individual does not survive
            if (uniform(rng_r) < mu(d, v))
            {
                metapop_iter->breeders[breeder_idx] = metapop_iter->breeders.back();

                metapop_iter->breeders.pop_back();

                // reduce index by 1
                // to now evaluate the survival of the individual at the back
                // of the stack that was copied to its new position at the 
                // middle of the stack.
                --breeder_idx;
            }
        }
    }
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

    // auxiliary variable to calculate group level fecundity
    double group_level_fecundity;

    // calculate a mean fecundity distribution
    for (auto metapop_iter = metapopulation.begin();
            metapop_iter != metapopulation.end();
            ++metapop_iter)
    {
        group_level_fecundity = 0.0;

        // calculate fecundity for each group
        // dependent on individual vigilance values
        for (auto breeder_iter = metapop_iter->breeders.begin();
                breeder_iter != metapop_iter->breeders.end();
                ++breeder_iter)
        {

            // calculate 1 - v^x
            group_level_fecundity += 1.0 - std::pow(breeder_iter->v[0] + breeder_iter->v[1], param.fecundity_power);
        }

        // add this value to the list of fecundities
        group_level_fecundities.push_back(group_level_fecundity);
    }

    // make a probability distribution
    std::discrete_distribution<unsigned> group_level_fecundity_distribution(
            group_level_fecundities.begin(), 
            group_level_fecundities.end());

    // TODO: let's park the use of the distribution for a while
    // which patch are we going to use? 
    unsigned patch_idx = group_level_fecundity_distribution(rng_r);


} // end StressSocial::reproduce()

// write parameters to file
void StressSocial::write_parameters(std::ofstream &data_file) 
{
    data_file << std::endl
        << std::endl
        << "seed;" << seed << ";" << std::endl
        << "time_step;" << time_step << ";" << std::endl
        << "dispersal;" << param.d << ";" << std::endl
        << "npatches;" << param.npatches << ";" << std::endl
        << "s;" << param.s << ";" << std::endl
        << "p_attack;" << param.p_attack << ";" << std::endl
        << "fecundity_power;" << param.fecundity_power << ";" << std::endl
        << "hmin;" << param.hmin << ";" << std::endl
        << "hmax;" << param.hmax << ";" << std::endl
        << "dmax;" << param.dmax << ";" << std::endl
        << "survival_power;" << param.survival_power << ";" << std::endl
        << "init_v;" << param.init_v << ";" << std::endl
        << "init_stress_hormone_level;" << param.init_stress_hormone_level << ";" << std::endl
        << "mu_baseline;" << param.mu_baseline << ";" << std::endl
        << "mu_stress_influx;" << param.mu_stress_influx << ";" << std::endl
        << "mu_vigilance_influx;" << param.mu_vigilance_influx << ";" << std::endl
        << "mu_removal;" << param.mu_removal << ";" << std::endl
        << "mu_v;" << param.mu_v << ";" << std::endl
        << "file_name;" << param.file_name << ";" << std::endl;
}

// update the stress hormone level for each individual
void StressSocial::update_stress_hormone()
{
    double stress_hormone_tplus1,r, stress_influx, vigilance_influx, baseline_influx;
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
            baseline_influx = breeder_iter->baseline_influx[0] + breeder_iter->baseline_influx[1];
            stress_influx = breeder_iter->stress_influx[0] + breeder_iter->stress_influx[1];
            vigilance_influx = breeder_iter->vigilance_influx[0] + breeder_iter->vigilance_influx[1];

            stress_hormone_tplus1 = (1.0 - r) * breeder_iter->stress_hormone + 
                baseline_influx +
                breeder_iter->is_attacked * stress_influx + 
                vigilance_influx * metapop_iter->V;

            // undo the is_attacked variable, ready for the next time step
            breeder_iter->is_attacked = false;
        }
    }
} // update_stress_hormone()

// updates damage each time step
void StressSocial::update_damage()
{

}

