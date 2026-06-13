#ifndef PDESOLVER_SOLVER_DRIVER_TRANSIENT_HPP
#define PDESOLVER_SOLVER_DRIVER_TRANSIENT_HPP

#include "solver/stage/Stage.hpp"
#include "solver/timestepper/TimeStepper.hpp"

namespace pdesolver {
	namespace solver {
		namespace driver {

			template<stage::Stage StageType, typename TimeStepperType>
			class Transient {
			public:

				bool solve(StageType& stage, TimeStepperType& stepper) {

					bool converged = true;

					while (!stepper.finished()) {
						
						stage.initialize();
						stage.assemble();
						converged = stage.solve();
						stepper.advance();
						stage.finalize();

					}

					return converged;
				}

			}; // class Transient

		} // namespace driver
	} // namespace solver
} // namespace pdesolver

#endif
