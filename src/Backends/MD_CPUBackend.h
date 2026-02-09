/**
 * @file    MD_CPUBackend.h
 * @date    03/set/2010
 * @author  lorenzo
 *
 *
 */

#ifndef MD_CPUBACKEND_H_
#define MD_CPUBACKEND_H_

#include "MDBackend.h"
#include "MCMoves/VolumeMove.h"

#include <iomanip>

class BaseThermostat;

/**
 * @brief Manages a MD simulation on CPU. It supports NVE and NVT simulations
 */

class MD_CPUBackend: public MDBackend {
protected:
	std::shared_ptr<BaseThermostat> _thermostat;
	MovePtr _V_move;
	int _stress_tensor_avg_every;
	int _stress_tensor_counter;

	// thermostat introduced in https://journals.aps.org/pre/abstract/10.1103/PhysRevE.75.056707
	bool _use_builtin_langevin_thermostat = false;
	number _langevin_c1 = 0.;
	number _langevin_c2 = 0.;

	void _first_step();
	void _compute_forces();
	void _second_step();

	void _update_backend_info();


	void _sync_crooks_data();
	std::vector<BaseForce*> _crooks_forces_registry;

    /**
     * @brief Helper to handle file I/O for Crooks forces on CPU.
     */
    template <typename T>
    void _process_crooks_sync_cpu(T* cpu_force, llint current) {
        const int buffer_size = 100000;
        
        if (current % buffer_size == 1) {
            cpu_force->saved_last_step = false;
            return;
        }

        if (cpu_force && !cpu_force->saved_last_step) {
            std::ofstream outputFile(cpu_force->_file_path, std::ios::app);
            if (outputFile.is_open()) {
                number _running_force = 0., _running_ext = 0.;
                
                for (int i = 0; i < buffer_size; i++) {
                    _running_force += cpu_force->_force_buffer[i];
                    _running_ext += cpu_force->_extension_buffer[i];
                    
                    if ((i + 1) % (cpu_force->_sum_steps) == 0) {
                        outputFile << std::setprecision(12) 
                                   << _running_force / cpu_force->_sum_steps << " " 
                                   << _running_ext / cpu_force->_sum_steps << " " 
                                   << cpu_force->_sum_steps << std::endl;
                        _running_force = 0.; _running_ext = 0.;
                    } 
                }
                outputFile.close();
                cpu_force->saved_last_step = true; 
            }
        }
    }

public:
	MD_CPUBackend();
	virtual ~MD_CPUBackend();

	void init();
	void get_settings(input_file &inp);
	void sim_step();
	void activate_thermostat();
};

#endif /* MD_CPUBACKEND_H_ */
