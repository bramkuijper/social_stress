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
    // 1. all patches that are of type P need to have a visit by a predator
    // 2. predator samples x individuals to attack 
    // 3. predator attacks them, so this changes an individuals' state
    for (auto metapop_iter = metapopulation.begin();
            metapop_iter != metapopulation.end();
            ++metapop_iter)
    {
        if (metapop_iter->predator_patch)
        {
            // option 1: sample total number of individuals that will be attacked
            // say, 3 individuals are being attacked. 
            for (auto female_iter = metapop_iter->breeders[female].begin();
                    female_iter < metapop_iter->breeders[female].end();
                    ++female_iter)
            {
                female_iter->is_attacked = uniform(rng_r) < param.p_attack;
                // ther will be code here telling about 
                // changes in hormone level or whatever
            }
            
            for (auto male_iter = metapop_iter->breeders[male].begin();
                    male_iter < metapop_iter->breeders[male].end();
                    ++male_iter)
            {
                male_iter->is_attacked = uniform(rng_r) < param.p_attack;
                // ther will be code here telling about 
                // changes in hormone level or whatever
            }
        }
    }
} // end predator_visit()


// probability of surviving an attack given hormone level h
double StressSocial::attack_survival(double const h)
{
    return(a_h * std::pow(h / hormone_max), b_h);
}

void StressSocial::survival()
{
    // aux variable containing the total mortality prob for an individual
    double survival_prob;

    for (auto metapop_iter = metapopulation.begin();
            metapop_iter != metapopulation.end();
            ++metapop_iter)
    {
        // first clear the list of previous survivors
        metapop_iter->breeders_surviving[female].clear();
        metapop_iter->breeders_surviving[male].clear();

         // have each and every individual survive dependent on 
        // (i) whether it was attacked
        // (ii) its baseline mortality and 
        // (iii) the level of damage it has accumulated
        //
        for (auto female_iter = metapop_iter->breeders[female].begin();
                female_iter < metapop_iter->breeders[female].end();
                ++female_iter)
        {

             survival_prob = 1.0;

            if (female_iter.is_attacked)
            {
                survival_prob *= attack_survival(
                        female_iter.stress_hormone[0] + 
                        female_iter.stress_hormone[1]
                        );
            }

            survival_prob *= param.max_survival;
       
            // TODO think about damage

            // TODO this is where we stopped on March 7 2024 14:29
            if (uniform(rng_r) < survival_prob)
            {
                metapop_iter->breeders_surviving[female].push_back(female_iter);
            }
        }
    }
} // end survival()

