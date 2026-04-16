#ifndef PDESOLVER_SOLVERREPORT_HPP
#define PDESOLVER_SOLVERREPORT_HPP

#include "core/Types.hpp"

namespace pdesolver {
	namespace linalg {
		namespace solver {
			
			template<typename DataType>
			struct SolverReport {
				
				bool converged = false;
				Index iterations = 0;
				DataType initialResidual = 0.0;
				DataType finalResidual = 0.0;

				dataType reduction() const {
					if (intialResidual == 0.0) return 0.0;
					return finalResidual / initialResidual;
				}

			}; // struct SolverReport

		} // namespace solver
	} // namespace linalg
} // namespace pdesolver
