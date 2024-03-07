#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

class Individual
{
    private:

    public:
        bool is_attacked{false};
        double stress_hormone[2]{0.0,0.0};
        double damage{0.0};

        Individual(double const init_stress_hormone_level);


};

#endif
