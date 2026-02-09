/*
 * MovingCrooksTrap.cpp
 *
 *  Created on: 18/jul/2025
 *      Author: Flavio 
 */

#include "MovingCrooksTrap.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include <iomanip>

#include <fstream>   // For std::ofstream


MovingCrooksTrap::MovingCrooksTrap() :
				BaseForce() {

}

std::tuple<std::vector<int>, std::string> MovingCrooksTrap::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "rate", &_rate, 1);

	getInputString(&inp, "file_path", _file_path, 1);
    _sum_steps = 1; //default stiff_rate is 0
    getInputInt(&inp, "sum_steps", &_sum_steps, 0);

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

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "MovingCrooksTrap");
	std::string description = Utils::sformat("MovingCrooksTrap (stiff=%g, rate=%g, dir=%g,%g,%g, pos0=%g,%g,%g", _stiff, _rate, _direction.x, _direction.y, _direction.z, _pos0.x, _pos0.y, _pos0.z);

	return std::make_tuple(particle_ids, description);
}

LR_vector MovingCrooksTrap::value(llint step, LR_vector &pos) {
	LR_vector postrap;
	number x, y, z;

	postrap.x = _pos0.x + (_rate * step) * _direction.x;
	postrap.y = _pos0.y + (_rate * step) * _direction.y;
	postrap.z = _pos0.z + (_rate * step) * _direction.z;

	x = -_stiff * (pos.x - postrap.x);
	y = -_stiff * (pos.y - postrap.y);
	z = -_stiff * (pos.z - postrap.z);

	LR_vector val = LR_vector(x, y, z);

    _force_buffer[step%100000] = val * _direction/_direction.module();
    _extension_buffer[step%100000] = (_rate * step);

    return val;
}

number MovingCrooksTrap::potential(llint step, LR_vector &pos) {
	LR_vector postrap;

	postrap.x = _pos0.x + (_rate * step) * _direction.x;
	postrap.y = _pos0.y + (_rate * step) * _direction.y;
	postrap.z = _pos0.z + (_rate * step) * _direction.z;

	return (number) (0.5 * _stiff * (pos - postrap).norm());
}
