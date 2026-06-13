#ifndef PDESOLVER_SOLVER_DRIVER_STEADY_HPP
#define PDESOLVER_SOLVER_DRIVER_STEADY_HPP

#include "solver/stage/Stage.hpp"

namespace pdesolver {
	namespace solver {
		namespace driver {

			template<stage::Stage StageType>
			class Steady {
			public:

				bool solve(StageType& stage) {

					stage.initialize();
					stage.assemble();
					bool converged = stage.solve();
					stage.finalize();

					return converged;
				}

			}; // class Steady

		} // namespace driver
	} // namespace solver
} // namespace pdesolver

#endif
