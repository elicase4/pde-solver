#ifndef PDESOLVER_SOLVER_CONFIG_LINEARSOLVERCONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_LINEARSOLVERCONFIG_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace solver {
		namespace config {

			struct LinearSolverConfig {

				enum class Type {
					CG,
					GMRES,
					BiCGSTAB,
					LU
				};

				Type type = Type::CG;

				Real tolerance = 1e-10;

				Index maxIterations = 1000;
				
				bool matrixFree = false;

				Index krylovDim = 50;

			};

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif
