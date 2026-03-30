#ifndef LTCOMABSPOSTtrap_H_
#define LTCOMABSPOSTRAP_H_

#include "../BaseForce.h"
#include <vector>
#include <string>

/**
 * @brief Metadynamics bias potential acting on the absolute position (along a
 *        single Cartesian axis) of the centre of mass of a particle list.
 *
 * The collective variable (CV) is:
 *   cv = (1/N) * sum_i get_abs_pos(p_i)[axis]
 *
 * which is identical to what the COMPosition observable reports.
 *
 * Input parameters:
 *   particle_list   = <string>  comma-separated particle indices
 *   axis            = <string>  "x", "y", or "z"
 *   xmin            = <float>   left boundary of the potential grid (inclusive)
 *   xmax            = <float>   right boundary of the potential grid (inclusive)
 *   N_grid          = <int>     number of grid points
 *   potential_grid  = <string>  comma-separated potential values (N_grid values)
 *
 * Type string in force file: "meta_com_abs_pos_trap"
 */

class LTCOMAbsPosTrap : public BaseForce {
public:
	std::vector<int> _p_list;

	/// Map from axis string to component index (0,1,2)
	int _axis_index = 0;

	std::vector<BaseParticle *> _p_list_ptr;

	number xmin = 0.;
	number xmax = 10.;
	int N_grid = 0;
	number dX = 0.1;
	std::vector<number> potential_grid;

	LTCOMAbsPosTrap();
	virtual ~LTCOMAbsPosTrap() {}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	LR_vector value(llint step, LR_vector &pos) override;
	number potential(llint step, LR_vector &pos) override;
};

#endif /* LTCOMABSPOSTRAP_H_ */
