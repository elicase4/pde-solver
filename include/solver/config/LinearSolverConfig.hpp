#ifndef RESIDUUM_SOLVER_CONFIG_LINEARSOLVERCONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_LINEARSOLVERCONFIG_HPP

#include "core/Types.hpp"

namespace residuum {
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
} // namespace residuum

#endif
