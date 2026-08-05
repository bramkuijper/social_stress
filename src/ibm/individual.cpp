#include <algorithm>
#include "individual.hpp"


// main constructor
Individual::Individual(Parameters const &params) :
    v{
        params.vigilance ? params.init_v / 2.0 : 0.0,
        params.vigilance ? params.init_v / 2.0 : 0.0
    }
{
    // Initialise stress-independent baseline influx
    baseline_influx[0] = params.init_baseline_influx;
    baseline_influx[1] = params.init_baseline_influx;
    
    // Initialise evolvable starting hormone trait.
    hstart[0] = params.init_hstart;
    hstart[1] = params.init_hstart;
    
    // Initialise current hormone from hstart.
    stress_hormone =
        0.5 * (hstart[0] + hstart[1]);

    stress_hormone =
        std::clamp(
            stress_hormone,
            0.0,
            params.hmax
        );
    
    // Initialising removal alleles
    removal[0] = params.init_removal;
    removal[1] = params.init_removal;
    
    // Initialise the stress-response feedback trait.
    h1_S[0] = params.init_h1_S;
    h1_S[1] = params.init_h1_S;

    // Start outside an active post-stressor response.
    hx = 0.0;
    time_since_last_stressor = params.tmax_stress_influx;
}


// copy constructor
Individual::Individual(Individual const &other) :
    is_alive{other.is_alive},
    is_attacked{other.is_attacked},
    baseline_influx{other.baseline_influx[0],other.baseline_influx[1]},
    stress_influx{other.stress_influx[0], other.stress_influx[1]},
    vigilance_influx{other.vigilance_influx[0], other.vigilance_influx[1]},
    removal{other.removal[0],other.removal[1]},
    h1_S{other.h1_S[0], other.h1_S[1]},
    hstart{other.hstart[0], other.hstart[1]},
    hx{other.hx},
    time_since_last_stressor{other.time_since_last_stressor},   
    v{other.v[0],other.v[1]},
    damage{other.damage},
    stress_hormone{other.stress_hormone}
{}

// birth constructor
Individual::Individual(
                Individual const &mum,
                Individual const &dad,
                Parameters const &param,
                std::mt19937 &rng_r) 
{
    std::bernoulli_distribution segregator{0.5};
    
    // TABORSKY VALIDATION
    // Each allelic value is bound to the full hormone scale. Because expressed
    // phenotype is the mean of the two alleles, this allows the phenotype to span the full range [0, hmax]

    baseline_influx[0] = mutate(mum.baseline_influx[segregator(rng_r)], param.mu_baseline, param.sdmu, rng_r);
    baseline_influx[0] = std::clamp(baseline_influx[0], 0.0, param.hmax);

    baseline_influx[1] = mutate(dad.baseline_influx[segregator(rng_r)], param.mu_baseline, param.sdmu, rng_r);
    baseline_influx[1] = std::clamp(baseline_influx[1], 0.0, param.hmax);

    stress_influx[0] = mutate(mum.stress_influx[segregator(rng_r)], param.mu_stress_influx, param.sdmu, rng_r);
    stress_influx[0] = std::clamp(stress_influx[0], 0.0, param.hmax);

    stress_influx[1] = mutate(dad.stress_influx[segregator(rng_r)], param.mu_stress_influx, param.sdmu, rng_r);
    stress_influx[1] = std::clamp(stress_influx[1], 0.0, param.hmax);
    
    // TABORSKY VALIDATION:
    // Inherit and mutate h1_S, which controls negative feedback
    // on stress-induced hormone production.
    h1_S[0] = mutate(
        mum.h1_S[segregator(rng_r)],
        param.mu_h1_S,
        param.sdmu,
        rng_r
    );
    h1_S[0] = std::clamp(h1_S[0], 0.0, 1.0);
    
    h1_S[1] = mutate(
        dad.h1_S[segregator(rng_r)],
        param.mu_h1_S,
        param.sdmu,
        rng_r
    );
    h1_S[1] = std::clamp(h1_S[1], 0.0, 1.0);
    
    
    if (param.vigilance)
    {
        vigilance_influx[0] = mutate(mum.vigilance_influx[segregator(rng_r)], param.mu_vigilance_influx, param.sdmu, rng_r);
        vigilance_influx[0] = std::clamp(vigilance_influx[0], 0.0, param.hmax);
      
        vigilance_influx[1] = mutate(dad.vigilance_influx[segregator(rng_r)], param.mu_vigilance_influx, param.sdmu, rng_r);
        vigilance_influx[1] = std::clamp(vigilance_influx[1], 0.0, param.hmax);
    
    }
      else
      {
      
        vigilance_influx[0] = 0.0;
        vigilance_influx[1] = 0.0;
      
      }
      
    removal[0] = mutate(mum.removal[segregator(rng_r)], param.mu_removal, param.sdmu, rng_r);
    removal[0] = std::clamp(removal[0],param.min_removal,1.0);

    removal[1] = mutate(dad.removal[segregator(rng_r)], param.mu_removal, param.sdmu, rng_r);
    removal[1] = std::clamp(removal[1],param.min_removal,1.0);
    
    // Inherit and mutate the starting hormone trait.
    hstart[0] = mutate(
        mum.hstart[segregator(rng_r)],
        param.mu_hstart,
        param.sdmu,
        rng_r
    );
    hstart[0] =
        std::clamp(hstart[0], 0.0, param.hmax);
    
    hstart[1] = mutate(
        dad.hstart[segregator(rng_r)],
        param.mu_hstart,
        param.sdmu,
        rng_r
    );
    hstart[1] =
        std::clamp(hstart[1], 0.0, param.hmax);

    // If vigilance disabled, vigilance alleles fixed at zero
    if (param.vigilance)
    {
    
      v[0] = mutate(mum.v[segregator(rng_r)], param.mu_v, param.sdmu, rng_r);
      v[0] = std::clamp(v[0], 0.0, 1.0);
      v[1] = mutate(dad.v[segregator(rng_r)], param.mu_v, param.sdmu, rng_r);
      v[1] = std::clamp(v[1], 0.0, 1.0);
    }
    else
    {
      v[0] = 0.0;
      v[1] = 0.0;
    }    
    
      // EG NOTE: Evolvable baseline vigilance alleles (a_v).
      // No b_stress?vigilance trait yet – vigilance is baseline-only.
      // Expressed vigilance is calculated via effective_vigilance()
      // = 0.5*(v[0] + v[1]) clamped to [0,1] in stress_social.cpp.

      // TABORSKY VALIDATION:
      // New offspring begin at their inherited hstart phenotype,
      // matching the Taborsky stress model.
      stress_hormone =
          0.5 * (hstart[0] + hstart[1]);
      
      stress_hormone =
          std::clamp(
              stress_hormone,
              0.0,
              param.hmax
          );

    // New offspring start outside an active post-stressor response.
    // These state variables will be updated when an attack occurs.
    hx = 0.0;
    time_since_last_stressor = param.tmax_stress_influx;


    // damage is 0 as per the default
    

} // birth constructor
  
void Individual::operator=(Individual const &other)
{
    is_alive = other.is_alive;
    is_attacked = other.is_attacked;
    damage = other.damage;
    stress_hormone = other.stress_hormone;
    hx = other.hx;
    time_since_last_stressor = other.time_since_last_stressor;

    for (unsigned allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        baseline_influx[allele_idx] = other.baseline_influx[allele_idx];
        stress_influx[allele_idx] = other.stress_influx[allele_idx];
        vigilance_influx[allele_idx] = other.vigilance_influx[allele_idx];
        removal[allele_idx] = other.removal[allele_idx];
        v[allele_idx] = other.v[allele_idx];
        h1_S[allele_idx] = other.h1_S[allele_idx];
        hstart[allele_idx] = other.hstart[allele_idx];
    }
} // end operator=()

// mutate according to a continuum-of-alleles model
double Individual::mutate(
        double to_mutate,
        double const mutation_prob,
        double const mutation_sd,
        std::mt19937 &rng_r)
{
    std::uniform_real_distribution<double> uniform{0.0,1.0};

    if (uniform(rng_r) < mutation_prob)
    {
        std::normal_distribution<double> mutational_effect_size{0.0,mutation_sd};

        // our new allelic value
        to_mutate += mutational_effect_size(rng_r);
    }

    return(to_mutate);
} // end mutate 


