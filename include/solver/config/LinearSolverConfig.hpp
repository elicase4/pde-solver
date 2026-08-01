#ifndef PDESOLVER_SOLVER_CONFIG_LINEARSOLVERCONFIG_HPP
#define PDESOLVER_SOLVER_CONFIG_LINEARSOLVERCONFIG_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace solver {
		namespace config {

			struct PreconditionerConfig {

				enum class Type {
					Identity
				};

				Type type = Type::Identity;

			}; // struct PreconditionerConfig

			struct LinearSolverConfig {

				enum class Type {
					CG,
					GMRES,
					BiCGSTAB,
					LU
				};

				enum class OperatorType {
					CSR,
					FEM
				};

				Type type = Type::CG;

				OperatorType operatorType = OperatorType::CSR;

				PreconditionerConfig preconditioner;

				Real tolerance = 1e-10;

				Index maxIterations = 1000;

				Index krylovDim = 50;

			};

		} // namespace config
	} // namespace solver
} // namespace pdesolver

#endif
