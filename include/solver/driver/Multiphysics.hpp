#ifndef PDESOLVER_SOLVER_DRIVER_MULTIPHYSICS_HPP
#define PDESOLVER_SOLVER_DRIVER_MULTIPHYSICS_HPP

#include <tuple>
#include <cmath>

#include "core/Types.hpp"
#include "solver/stage/Stage.hpp"

namespace pdesolver {
	namespace solver {
		namespace driver {

			template<typename ResidualFunctor, stage::Stage... StageTypes>
			class MultiphysicsDriver {
			public:

				struct Config {
					Index maxOuterIterations = 50;
					Real absoluteTolerance = 1e-8;
					Real relativeTolerance = 1e-6;
				}; // struct Config
	
				explicit MultiphysicsDriver(Config cfg, ResidualFunctor residual, StageTypes&... stages) : cfg_(cfg), residual_(residual), stages_(stages) {}

				bool solve() {

					Real residual0 = -1.0;

					for (Index iter = 0; iter < cfg_.maxOuterIteration; ++iter) {

						// segregated solve for all stages in specified order
						bool allConverged = sweepAll(std::index_sequence_for<StageTypes...>{});
						if (!allConverged) {
							return false;
						}

						// evaluate inter-field coupling residual
						Real res = residual_(iter);
						if (iter == 0) {
							residual0 = res;
						}

						// check absolute and relative convergence
						if (res < cfg_.absoluteTolerance) {
							return true;
						}

						if ((residual0 > 0.0 && (res / residual0)) < cfg_.relativeTolerance) {
							return true;
						}

					}

					return false;

				}

			private:

				template<Index... Is>
				bool sweepAll(std::index_sequence<Is...>) {
					return (... && solveStage(std::get<Is>(stages_)) );
				}

				template<typename S>
				bool solveStage(S& stage) {
					stage.assemble();
					return stage.solve();
				}

				Config cfg_;
				ResidualFunctor residual_;
				std::tuple<StageType&...> stages_;

			}; // class MultiphysicsDriver

		} // namespace driver
	} // namespace solver
} // namespace pdesolver

#endif
