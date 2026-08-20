#ifndef PDESOLVER_UTILS_LOGGING_SOLVER_LOGGER_HPP
#define PDESOLVER_UTILS_LOGGING_SOLVER_LOGGER_HPP

#include <utility>
#include <variant>
#include <vector>

#include "core/Types.hpp"
#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/solver/ConsoleLogger.hpp"

namespace pdesolver {
	namespace utils {
		namespace logging {

			// Runtime choice between ConsoleLogger and NullLogger for the per-solve
			// iteration/residual channel -- mirrors driver::Logger's variant pattern so
			// LinearSolverFactory can honor LoggingConfig.solver.type.
			class SolverLogger {
			public:

				using Variant = std::variant<ConsoleLogger, NullLogger>;

				explicit SolverLogger(Variant impl) : impl_(std::move(impl)) {}

				template<typename DataType>
				void log(Index iter, const std::vector<DataType>& perDOFAbs, DataType flopsThisIter = DataType(0)) const {
					std::visit([&](const auto& l){ l.log(iter, perDOFAbs, flopsThisIter); }, impl_);
				}

				template<typename DataType>
				std::vector<DataType> computePerDOFNorms(const DataType* r, Index totalSize) const {
					return std::visit([&](const auto& l){ return l.template computePerDOFNorms<DataType>(r, totalSize); }, impl_);
				}

				void summary(bool converged) const {
					std::visit([&](const auto& l){ l.summary(converged); }, impl_);
				}

				template<typename Args>
				void event(Args&& msg) const {
					std::visit([&](const auto& l){ l.event(std::forward<Args>(msg)); }, impl_);
				}

			private:
				Variant impl_;

			}; // class SolverLogger

		} // namespace logging
	} // namespace utils
} // namespace pdesolver

#endif
