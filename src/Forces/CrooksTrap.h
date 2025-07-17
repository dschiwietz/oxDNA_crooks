/**
* @file    CrooksTrap.h
* @date    25/Jun/2025
* @author  Dominik
*
*/

#ifndef CROOKSTRAP_H_
#define CROOKSTRAP_H_

#include "BaseForce.h"

using namespace std;

class CrooksTrap: public BaseForce {
private:
    int _particle;
    int _ref_id;
    static number _work_buffer[100000];
    static bool saved_last_step;
    llint last_step;

public:
    BaseParticle * _p_ptr;
    number _r0;
    number _rate;
    number _stiff_rate;
    string _file_path;
    int _sum_steps;
    bool PBC;

    CrooksTrap();
    virtual ~CrooksTrap() {
    }

    std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

    virtual LR_vector value(llint step, LR_vector &pos);
    virtual number potential(llint step, LR_vector &pos);

protected:
    LR_vector _distance(LR_vector u, LR_vector v);
};

#endif // CROOKSTRAP_H
 