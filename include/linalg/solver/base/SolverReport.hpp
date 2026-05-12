#ifndef PDESOLVER_SOLVERREPORT_HPP
#define PDESOLVER_SOLVERREPORT_HPP

#include <vector>

#include "core/Types.hpp"

namespace pdesolver {
	namespace linalg {
		namespace solver {
			
			template<typename VectorType>
			struct SolverReport {
				
				using DataType = typename VectorType::value_type;

				bool converged = false;
				Index iterations = 0;
				DataType initialResidual = 0.0; // || r_0 ||
				DataType finalResidual = 0.0; // || r_k ||
				DataType finalResidualRel = 0.0; // finalResidual / initialResidual
				std::vector<DataType> perFieldResidual;

				DataType reduction() const {
					return (initialResidual > DataType(0)) ? (finalResidual / initialResidual) : DataType(0);
				}

			}; // struct SolverReport

		} // namespace solver
	} // namespace linalg
} // namespace pdesolver

#endif
