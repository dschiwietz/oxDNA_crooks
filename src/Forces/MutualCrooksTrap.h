/**
* @file    MutualCrooksTrap.h
* @date    25/Jun/2025
* @author  Dominik
*
*/

#ifndef MutualCrooksTrap_H_
#define MutualCrooksTrap_H_

#include "BaseForce.h"
#include "MutualCrooksTrap.h"

using namespace std;

class MutualCrooksTrap: public BaseForce {
private:
    int _particle;
    int _ref_id;

public:
    BaseParticle * _p_ptr;
    number _r0;
    number _rate;
    number _stiff_rate;
    string _file_path;
    number _force_buffer[100000];
    number _extension_buffer[100000];
    bool saved_last_step;
    llint last_step;
    int _sum_steps;
    bool PBC;

    MutualCrooksTrap();
    virtual ~MutualCrooksTrap() {
    }

    std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

    virtual LR_vector value(llint step, LR_vector &pos);
    virtual number potential(llint step, LR_vector &pos);

protected:
    LR_vector _distance(LR_vector u, LR_vector v);
};

#endif // MutualCrooksTrap_H
 