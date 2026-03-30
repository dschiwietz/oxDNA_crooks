/*
 * LTCOMAbsPosTrap.cpp
 *
 * Metadynamics bias force whose collective variable (CV) is the absolute
 * position of the centre of mass of a particle list along one Cartesian axis.
 *
 * This mirrors what the COMPosition observable computes, making the CV
 * directly comparable to that observable's output.
 *
 * Author: dominik
 */

#include "LTCOMAbsPosTrap.h"

#include "meta_utils.h"
#include "../../Particles/BaseParticle.h"

#include <string>
#include <iostream>
#include <cmath>

// ── Helpers (file-local, same pattern as LTCOMTrap) ──────────────────────────

namespace {

/// Linear interpolation of the potential at position x.
inline number interpolate_potential_1d(const number x,
                                       const number dX,
                                       const number xmin,
                                       const std::vector<number> &grid) {
	number ix_left_f = std::floor((x - xmin) / dX);
	int ix_left  = static_cast<int>(ix_left_f);
	int ix_right = ix_left + 1;

	number x_left  = xmin + ix_left * dX;
	number x_right = x_left + dX;

	return (x_right - x) / dX * grid[ix_left] + (x - x_left) / dX * grid[ix_right];
}

/// Finite-difference derivative → force along the CV axis.
inline number get_axis_force_1d(const number x,
                                 const number dX,
                                 const number xmin,
                                 const std::vector<number> &grid) {
	int ix_left  = static_cast<int>(std::floor((x - xmin) / dX));
	int ix_right = ix_left + 1;
	// F = -dV/dx, first-order finite difference
	return -(grid[ix_right] - grid[ix_left]) / dX;
}

} // anonymous namespace

// ── LTCOMAbsPosTrap ──────────────────────────────────────────────────────────

LTCOMAbsPosTrap::LTCOMAbsPosTrap() : BaseForce() {
	xmin   = 0.;
	xmax   = 10.;
	N_grid = 0;
	dX     = 0.1;
}

std::tuple<std::vector<int>, std::string> LTCOMAbsPosTrap::init(input_file &inp) {
	BaseForce::init(inp);

	// Particle list
	std::tie(_p_list, _p_list_ptr) = meta::get_particle_lists(
		inp, "particle_list", CONFIG_INFO->particles(), "LTCOMAbsPosTrap particle_list");

	// Axis
	std::string axis_str;
	getInputString(&inp, "axis", axis_str, 1);
	if(axis_str == "x") {
		_axis_index = 0;
	} else if(axis_str == "y") {
		_axis_index = 1;
	} else if(axis_str == "z") {
		_axis_index = 2;
	} else {
		throw oxDNAException("LTCOMAbsPosTrap: invalid axis '%s' (must be 'x', 'y', or 'z')", axis_str.c_str());
	}

	// Grid parameters
	getInputNumber(&inp, "xmin",   &xmin,   1);
	getInputNumber(&inp, "xmax",   &xmax,   1);
	getInputInt   (&inp, "N_grid", &N_grid, 1);
	potential_grid.reserve(N_grid);

	std::string potential_string;
	getInputString(&inp, "potential_grid", potential_string, 1);
	potential_grid = meta::split_to_numbers(potential_string, ",");

	if((int) potential_grid.size() != N_grid) {
		throw oxDNAException(
			"LTCOMAbsPosTrap: potential_grid has %zu values but N_grid = %d",
			potential_grid.size(), N_grid);
	}

	dX = (xmax - xmin) / (N_grid - 1.0);

	std::string description = Utils::sformat(
		"LTCOMAbsPosTrap (axis=%s, xmin=%g, xmax=%g, N_grid=%d)",
		axis_str.c_str(), xmin, xmax, N_grid);

	return std::make_tuple(_p_list, description);
}

LR_vector LTCOMAbsPosTrap::value(llint step, LR_vector &pos) {
	// Compute absolute COM position along the chosen axis
	number com_axis = 0.;
	for(auto p : _p_list_ptr) {
		com_axis += CONFIG_INFO->box->get_abs_pos(p)[_axis_index];
	}
	com_axis /= static_cast<number>(_p_list_ptr.size());

	int ix_left  = static_cast<int>(std::floor((com_axis - xmin) / dX));
	int ix_right = ix_left + 1;

	LR_vector force(0., 0., 0.);

	if(ix_left < 0 || ix_right > N_grid - 1) {
		std::cout << "LTCOMAbsPosTrap: off grid! com_axis=" << com_axis << std::endl;
	} else {
		number F_axis = get_axis_force_1d(com_axis, dX, xmin, potential_grid);
		// Distribute equally over all particles in the COM list
		number F_per_particle = F_axis / static_cast<number>(_p_list_ptr.size());

		if(_axis_index == 0)      force = LR_vector(F_per_particle, 0., 0.);
		else if(_axis_index == 1) force = LR_vector(0., F_per_particle, 0.);
		else                      force = LR_vector(0., 0., F_per_particle);
	}

	return force;
}

number LTCOMAbsPosTrap::potential(llint step, LR_vector &pos) {
	// Compute absolute COM position along the chosen axis
	number com_axis = 0.;
	for(auto p : _p_list_ptr) {
		com_axis += CONFIG_INFO->box->get_abs_pos(p)[_axis_index];
	}
	com_axis /= static_cast<number>(_p_list_ptr.size());

	int ix_left  = static_cast<int>(std::floor((com_axis - xmin) / dX));
	int ix_right = ix_left + 1;

	if(ix_left < 0 || ix_right > N_grid - 1) {
		std::cout << "LTCOMAbsPosTrap: off grid! com_axis=" << com_axis << std::endl;
		return 0.;
	}

	number V = interpolate_potential_1d(com_axis, dX, xmin, potential_grid);
	// Potential is per particle (consistent with force distribution)
	return V / static_cast<number>(_p_list_ptr.size());
}
