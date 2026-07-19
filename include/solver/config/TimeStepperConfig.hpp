#ifndef PDESOLVER_SOLVER_CONFIG_TIMESTEPPERCONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_TIMESTEPPERCONFIG_HPP

#include "core/Types.hpp"
#include "solver/config/NonlinearSolverConfig.hpp"
#include "solver/config/LinearSolverConfig.hpp"

namespace pdesolver {
	namespace solver {
		namespace config {
			
			struct TimeStepperConfig {

				enum class Type {
					ForwardEuler,
					BackwardEuler,
					GeneralizedAlpha,
					RK4
				}; // enum class Type
	
				Type type = Type::ForwardEuler;
				
				// Time interval
				Real t0 = 0.0;
				Real tf = 1.0;
				Real dt = 1e-3;

				// Generalized alpha spectral radius
				// Real rhoInf = 1 // no numerical dissipation
				// Real rhoInf = 0 // maximum high-frequency dissipation
				Real rhoInf = 0.5;

				// Solver config
				NonlinearSolverConfig nonlinearSolver;
				LinearSolverConfig linearSolver;

			}; // struct TimeStepperConfig

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif
