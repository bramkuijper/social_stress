// heart of the stress social code

#include "stress_social.hpp"

// constructor function
StressSocial::StressSocial(Parameters const &parvals) :
    param{parvals}
    ,rd{} // initialize random device
    ,seed{rd()} // initialize seed
    ,rng_r{seed} // initialize the random number generator
    ,uniform{0.0,1.0} // initialize the uniform dist between 0 and 1
    ,patch_sampler{0, param.npatches - 1} // initialize uniform distribution to sample patch indices from
    ,metapopulation(param.npatches, Patch()) // initialize the metapopulation
    ,data_file{param.file_name}
{
    // make some patches P and some NP
    initialize_patches();


    // now run the bloody thing
    for (time_step = 0; time_step < param.max_time; ++time_step)
    {
        // do something!
        // effectively, we want the predator to visit some patches
        // and to attack individuals there.
        predator_visit();
        survival();
    }
}

// go over all the patches and initialize them as type NP or P
void StressSocial::initialize_patches()
{
    // calculate probability of encountering a predator on a patch
    double prob_P = param.s[NP] / (param.s[NP] + param.s[P]);

    for (auto patch_iterator = metapopulation.begin();
            patch_iterator != metapopulation.end();
            ++patch_iterator)
    {
        patch_iterator->predator_patch = uniform(rng_r) < prob_P;
    }
}

void StressSocial::write_data() 
{
    data_file << time_step << std::endl;
}

void StressSocial::predator_visit()
{
    double V; // auxiliary variable reflecting 
              // whether at least a single individual is vigilant

    // 1. all patches that are of type P need to have a visit by a predator
    // 2. predator samples x individuals to attack 
    // 3. predator attacks them, so this changes an individuals' state
    for (auto metapop_iter = metapopulation.begin();
            metapop_iter != metapopulation.end();
            ++metapop_iter)
    {
        if (metapop_iter->predator_patch)
        {
            // check whether at least one individual is vigilant
            V = calculate_group_vigilance(metapop_iter);

            // calculate the probability that nobody is vigilant
            if (uniform(rng_r) < 1.0 - V)
            {
                // then sample which individual will die
                random_breeder_idx = take_random_breeder(rng_r);

                // overwrite breeder at position with final breeder in stack
                metapopulation.breeders[random_breeder_idx] = metapopulation.breeders.back();

                // delete final element
                metapopulation.breeders.pop_back();
            }
        }
    }
} // end predator_visit()

double StressSocial::calculate_group_vigilance(Patch const &current_patch)
{
    double prod_1minusv = 1.0;
    // go over all patches and calculate vigilance
    for (auto breeder_iter = current_patch.breeders.begin();
            breeder_iter != current_patch.breeders.end();
            ++breeder_iter)
    {
        prod_1minusv = prod_1minusv * (1.0 - breeder_iter->v);
    }

    // return 1 - (1-v)^n
    return 1.0 - prod_1minusv;
} // end calculate_group_vigilance()


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
        for (auto breeder_iter = metapop_iter.breeders.begin();
                breeder_iter != metapop_iter.breeders.end();
                ++breeder_iter)
        {

            // calculate 1 - v^x
            group_level_fecundity += 1.0 - std::pow(breeder_iter->v, par.fecundity_power)
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

