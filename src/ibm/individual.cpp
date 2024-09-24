#include <algorithm>
#include "individual.hpp"


// main constructor
Individual::Individual(Parameters const &params) :
    v{params.init_v/2, params.init_v/2},
    stress_hormone{params.init_stress_hormone_level}
{}

// copy constructor
Individual::Individual(Individual const &other) :
    is_alive{other.is_alive},
    is_attacked{other.is_attacked},
    baseline_influx{other.baseline_influx[0],other.baseline_influx[1]},
    stress_influx{other.stress_influx[0], other.stress_influx[1]},
    vigilance_influx{other.vigilance_influx[0], other.vigilance_influx[1]},
    removal{other.removal[0],other.removal[1]},
    v{other.v[0],other.v[1]},
    damage{other.damage},
    stress_hormone{other.stress_hormone}
{}

// birth constructor
Individual::Individual(Individual const &mum,
                Individual const &dad,
                Parameters const &param,
                std::mt19937 &rng_r)
{
    std::bernoulli_distribution segregator{0.5};

    baseline_influx[0] = mutate(mum.baseline_influx[segregator(rng_r)], param.mu_baseline, param.sdmu, rng_r);
    baseline_influx[0] = std::clamp(baseline_influx[0], 0.0, param.hmax/2.0);

    baseline_influx[1] = mutate(dad.baseline_influx[segregator(rng_r)], param.mu_baseline, param.sdmu, rng_r);
    baseline_influx[1] = std::clamp(baseline_influx[1], 0.0, param.hmax/2.0);

    stress_influx[0] = mutate(mum.stress_influx[segregator(rng_r)], param.mu_stress_influx, param.sdmu, rng_r);
    stress_influx[0] = std::clamp(stress_influx[0], 0.0, param.hmax/2.0);

    stress_influx[1] = mutate(dad.stress_influx[segregator(rng_r)], param.mu_stress_influx, param.sdmu, rng_r);
    stress_influx[1] = std::clamp(stress_influx[1], 0.0, param.hmax/2.0);
    
    vigilance_influx[0] = mutate(mum.vigilance_influx[segregator(rng_r)], param.mu_vigilance_influx, param.sdmu, rng_r);
    vigilance_influx[0] = std::clamp(vigilance_influx[0], 0.0, param.hmax/2.0);

    vigilance_influx[1] = mutate(dad.vigilance_influx[segregator(rng_r)], param.mu_vigilance_influx, param.sdmu, rng_r);
    vigilance_influx[1] = std::clamp(vigilance_influx[1], 0.0, param.hmax/2.0);
    
    removal[0] = mutate(mum.removal[segregator(rng_r)], param.mu_removal, param.sdmu, rng_r);
    removal[0] = std::clamp(removal[0], 0.0, 1.0);

    removal[1] = mutate(dad.removal[segregator(rng_r)], param.mu_removal, param.sdmu, rng_r);
    removal[1] = std::clamp(removal[1], 0.0, 1.0);

    v[0] = mutate(mum.v[segregator(rng_r)], param.mu_v, param.sdmu, rng_r);
    v[0] = std::clamp(v[0], 0.0, 1.0);
    v[1] = mutate(dad.v[segregator(rng_r)], param.mu_v, param.sdmu, rng_r);
    v[1] = std::clamp(v[1], 0.0, 1.0);


    // updating the state variables
    stress_hormone = (baseline_influx[0] + baseline_influx[1])/(removal[0] + removal[1]);

    // damage is 0 as per the default

} // birth constructor
  
void Individual::operator=(Individual const &other)
{
    is_alive = other.is_alive;
    is_attacked = other.is_attacked;
    damage = other.damage;
    stress_hormone = other.stress_hormone;

    for (unsigned allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        baseline_influx[allele_idx] = other.baseline_influx[allele_idx];
        stress_influx[allele_idx] = other.stress_influx[allele_idx];
        vigilance_influx[allele_idx] = other.vigilance_influx[allele_idx];
        removal[allele_idx] = other.removal[allele_idx];
        v[allele_idx] = other.v[allele_idx];
    }
} // end operator=()

// mutate according to a continuum-of-alleles model
double Individual::mutate(double to_mutate,
        double const mutation_prob,
        double const mutation_sd,
        std::mt19937 &rng_r)
{
    std::uniform_real_distribution<double> uniform{0.0,1.0};

    if (mutation_prob < uniform(rng_r))
    {
        std::normal_distribution<double> mutational_effect_size{0.0,mutation_sd};

        // our new allelic value
        to_mutate += mutational_effect_size(rng_r);
    }

    return(to_mutate);
} // end mutate 


