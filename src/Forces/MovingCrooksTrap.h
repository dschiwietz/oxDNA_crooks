/**
 * @file    MovingCrooksTrap.h
 * @date    18/oct/2014
 * @author  Flavio 
 *
 */

#ifndef MovingCrooksTrap_H_
#define MovingCrooksTrap_H_

#include "BaseForce.h"

using namespace std;

class MovingCrooksTrap: public BaseForce {
public:
	MovingCrooksTrap();
	virtual ~MovingCrooksTrap() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);

	string _file_path;
	number _force_buffer[100000];
    number _extension_buffer[100000];
    bool saved_last_step;
    llint last_step;
	int _sum_steps;
};

#endif // MovingCrooksTrap_H
