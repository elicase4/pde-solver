#ifndef PDESOLVER_SOLVER_DRIVER_TRANSIENT_HPP
#define PDESOLVER_SOLVER_DRIVER_TRANSIENT_HPP

#include "solver/stage/Stage.hpp"
#include "solver/timestepper/TimeStepperRunner.hpp"

namespace pdesolver {
	namespace solver {
		namespace driver {

			// Owns only the outer loop -- one-time initialize()/finalize() plus stepping until
			// the timestepper says it's done. Each step's assemble/solve/advance/output lives in
			// the timestepper (see TimeStepperRunner). Mirrors Driver::Steady.
			template<stage::Stage StageType>
			class Transient {
			public:

				bool solve(StageType& stage, timestepper::TimeStepperRunner& stepper) {

					stage.initialize();

					while (!stepper.finished()) {
						if (!stepper.step()) {
							stage.finalize();
							return false;
						}
					}

					stage.finalize();
					return true;

				}

			}; // class Transient

		} // namespace driver
	} // namespace solver
} // namespace pdesolver

#endif
