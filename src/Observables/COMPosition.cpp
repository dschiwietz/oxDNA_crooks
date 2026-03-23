/*
 * Distance.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#include "COMPosition.h"
#include "../Utilities/Utils.h"

COMPosition::COMPosition() {
	_PBC = true;
}

COMPosition::~COMPosition() {

}

void COMPosition::init() {
	BaseObservable::init();

	int N = _config_info->N();
	std::vector<BaseParticle *> &particles = _config_info->particles();

	std::vector<int> p1_indexes = Utils::get_particles_from_string(particles, _p1_string, "Distance observable");
	for(auto idx: p1_indexes) {
		_check_index(idx, N);
		_p1_list.insert(particles[idx]);
	}
	if (axis == "x"){
		index = 0;
	}else if (axis == "y"){
		index = 1;
	}else if (axis == "z"){
		index = 2;
	}else{
		throw oxDNAException("COMPosition: invalid axis %s", axis);
	}
	

}

std::string COMPosition::get_output_string(llint curr_step) {
	LR_vector dist;
	LR_vector p1_com, p2_com;
	for(auto particle : _p1_list) {
		p1_com += _config_info->box->get_abs_pos(particle);
	}
	p1_com /= _p1_list.size();

	return Utils::sformat("%14.4lf", p1_com[index]);
}

void COMPosition::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	// particle 1
	getInputString(&my_inp, "particle_1", _p1_string, 1);

	getInputBool(&my_inp, "PBC", &_PBC, 0);

	getInputString(&my_inp, "axis", axis, 1);
}
