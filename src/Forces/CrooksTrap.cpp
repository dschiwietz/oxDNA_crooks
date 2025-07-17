/*
* CrooksTrap.cpp
*
*  Created on: 25/Jun/2025
*      Author: Dominik 
*/

#include "CrooksTrap.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

#include <fstream>   // For std::ofstream

void appendBufferToFile(const std::string& filename, number* buffer, int step) {
    std::ofstream outputFile(filename, std::ios::app);

    if (outputFile.is_open()) {
        number _running = 0.;
        for (int i = 0; i < 100000; i++) {
            _running += buffer[i];
            if ((i+1) % (step) == 0){
                outputFile << _running << std::endl;
                _running = 0.;
            } 
        }
        outputFile.close();
    } else {
        std::cerr << "Error: Unable to open file '" << filename << "' for appending." << std::endl;
    }
}

number CrooksTrap::_work_buffer[100000] = {};
bool CrooksTrap::saved_last_step = false;

CrooksTrap::CrooksTrap() :
                BaseForce() {
    _ref_id = -2;
    _particle = -2;
    _p_ptr = NULL;
    _r0 = -1.;
    _rate = -1;
    _stiff_rate = -1;
    PBC = false;
}

std::tuple<std::vector<int>, std::string> CrooksTrap::init(input_file &inp) {
    BaseForce::init(inp);

    getInputInt(&inp, "particle", &_particle, 1);
    getInputInt(&inp, "ref_particle", &_ref_id, 1);
    getInputNumber(&inp, "r0", &_r0, 1);
    getInputNumber(&inp, "stiff", &_stiff, 1);
    getInputBool(&inp, "PBC", &PBC, 0);
    _rate = 0.f; //default rate is 0
    getInputNumber(&inp, "rate", &_rate, 0);
    _stiff_rate = 0.f; //default stiff_rate is 0
    getInputNumber(&inp, "stiff_rate", &_stiff_rate, 0);
    getInputString(&inp, "file_path", _file_path, 1);
    _sum_steps = 1; //default stiff_rate is 0
    getInputInt(&inp, "sum_steps", &_sum_steps, 0);

    int N = CONFIG_INFO->particles().size();
    if(_ref_id < 0 || _ref_id >= N) {
        throw oxDNAException("Invalid reference particle %d for Crooks Trap", _ref_id);
    }
    _p_ptr = CONFIG_INFO->particles()[_ref_id];

    if(_particle >= N || N < -1) {
        throw oxDNAException("Trying to add a CrooksTrap on non-existent particle %d. Aborting", _particle);
    }
    if(_particle == -1) {
        throw oxDNAException("Cannot apply CrooksTrap to all particles. Aborting");
    }

    std::string description = Utils::sformat("CrooksTrap (stiff=%g, stiff_rate=%g, r0=%g, rate=%g, ref_particle=%d, PBC=%d)", _stiff, _stiff_rate, _r0, _rate, _ref_id, PBC);

    return std::make_tuple(std::vector<int>{_particle}, description);
}

LR_vector CrooksTrap::_distance(LR_vector u, LR_vector v) {
    if(PBC) {
        return CONFIG_INFO->box->min_image(u, v);
    }
    else {
        return v - u;
    }
}

LR_vector CrooksTrap::value(llint step, LR_vector &pos) {
    LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr));
    LR_vector val = (dr / dr.module()) * (dr.module() - (_r0 + (_rate * step))) * (_stiff + (_stiff_rate * step));
    
    if (step > last_step) {
        if (step%100000 == 0 and step != 0 and !saved_last_step) {
            saved_last_step = true;
            appendBufferToFile(_file_path, _work_buffer, _sum_steps);
            for (int i = 0; i < 100000; ++i) {
                _work_buffer[i] = 0;
            }
        }else if (step%100000 == 1) {
            saved_last_step = false;
        }
        _work_buffer[step%100000] += val * dr/dr.module() * _rate;
        last_step = step;
    } 
    return val;
}

number CrooksTrap::potential(llint step, LR_vector &pos) {
    LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr));
    return pow(dr.module() - (_r0 + (_rate * step)), 2) * ((number) 0.5) * (_stiff + (_stiff_rate * step));
}
