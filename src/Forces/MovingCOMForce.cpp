/*
 * COMForce.cpp
 *
 *  Created on: 16 October 2025
 *      Author: dominik
 */

#include "MovingCOMForce.h"

#include "../Utilities/Utils.h"
#include "../Boxes/BaseBox.h"

#include <vector>

using namespace std;

MovingCOMForce::MovingCOMForce() {
	_rate = 0;
}

MovingCOMForce::~MovingCOMForce() {

}

std::tuple<std::vector<int>, std::string> MovingCOMForce::init(input_file &inp) {
	BaseForce::init(inp);

	getInputString(&inp, "com_list", _com_string, 1);
	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "rate", &_rate, 0);

	int tmpi;
	double tmpf[3];
	std::string strdir;
	getInputString(&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
	}
	_direction = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	_direction.normalize();

	getInputString(&inp, "pos0", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("Could not parse pos0 %s in external forces file. Aborting", strdir.c_str());
	}
	_pos0 = LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);


	auto com_indexes = Utils::get_particles_from_string(CONFIG_INFO->particles(), _com_string, "COMForce");
	for(auto it = com_indexes.begin(); it != com_indexes.end(); it++) {
		_com_list.insert(CONFIG_INFO->particles()[*it]);
	}

	std::string description = Utils::sformat("MovingCOMTrap (stiff=%g, rate=%g, dir=%g,%g,%g, pos0=%g,%g,%g)", _stiff, _rate, _direction.x, _direction.y, _direction.z, _pos0.x, _pos0.y, _pos0.z);
	return std::make_tuple(com_indexes, description);
}

void MovingCOMForce::_compute_coms(llint step) {
	if(step != _last_step) {
		_com = LR_vector(0, 0, 0);
		for(auto p : _com_list) {
			_com += CONFIG_INFO->box->get_abs_pos(p);
		}
		_com /= _com_list.size();

		_last_step = step;
	}
}

LR_vector MovingCOMForce::value(llint step, LR_vector &pos) {
	_compute_coms(step);
	LR_vector dist = (_pos0 + (_rate * step) * _direction - _com);
	number d_com = dist.module();
	number force = d_com * _stiff / _com_list.size();

	return dist * (force / d_com);
}

number MovingCOMForce::potential(llint step, LR_vector &pos) {
	_compute_coms(step);
	return 0.5 * _stiff * SQR((_pos0 + (_rate * step) * _direction - _com).module()) / _com_list.size();
}
