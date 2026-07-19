#ifndef PDESOLVER_SOLVER_DRIVER_TRANSIENT_HPP
#define PDESOLVER_SOLVER_DRIVER_TRANSIENT_HPP

#include "solver/stage/Stage.hpp"
#include "solver/timestepper/TimeStepper.hpp"

namespace pdesolver {
	namespace solver {
		namespace driver {

			template<stage::Stage StageType, typename TimeStepperType, typename VectorType>
			class Transient {
			public:

				bool solve(StageType& stage, TimeStepperType& stepper, VectorType& U, VectorType& U_prev) {

					while (!stepper.finished()) {
						
						stage.initialize();
						stage.assemble();

						if (!stage.solve()) {
							return false;
						}

						stepper.advance(U, U_prev);
						stage.finalize();

					}

					return true;
				}

			}; // class Transient

		} // namespace driver
	} // namespace solver
} // namespace pdesolver

#endif
