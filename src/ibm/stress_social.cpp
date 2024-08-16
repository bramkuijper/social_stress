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
        
        survive_damage();

        reproduce();
    }
} // end StressSocial constructor

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
    // 4. an individual can avoid attack dependent on its strrress response
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

                // we need to implement that individuals can flee the attack 
                // dependent on their stress hormone level h
                if (uniform(rng_r) < 1.0 - attack_survival(metapop_iter->breeders[random_breeder_idx].h))
                {
                    // overwrite breeder at position with final breeder in stack
                    metapop_iter->breeders[random_breeder_idx] = metapop_iter->breeders.back();

                    // delete final element
                    metapop_iter->breeders.pop_back();
                }
            }
        }
    }
} // end predator_visit()


// probability of surviving an attack given hormone level h
double StressSocial::attack_survival(double const h)
{
    return(pow(h/param.hmax, param.survival_power));
}

double StressSocial::calculate_group_vigilance(Patch const &current_patch)
{
    double prob_none_vigilant = 1.0;

    // go over all patches and calculate the total probability
    // that none of the individuals are vigilant
    for (auto breeder_iter = current_patch.breeders.begin();
            breeder_iter != current_patch.breeders.end();
            ++breeder_iter)
    {
        prob_none_vigilant = prob_none_vigilant * (1.0 - breeder_iter->v);
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
// TODO: we need to update damage levels each and every generation
void StressSocial::survive_damage_vigilance()
{
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
            d = metapop_iter->breeders[breeder_idx].d;
            v = metapop_iter->breeders[breeder_idx].v;

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

