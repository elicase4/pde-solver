#ifndef RESIDUUM_UTILS_LOGGING_SOLVER_LOGGER_HPP
#define RESIDUUM_UTILS_LOGGING_SOLVER_LOGGER_HPP

#include <utility>
#include <variant>
#include <vector>

#include "core/Types.hpp"
#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/solver/ConsoleLogger.hpp"

namespace residuum {
	namespace utils {
		namespace logging {

			class SolverLogger {
			public:

				using Variant = std::variant<ConsoleLogger, NullLogger>;

				explicit SolverLogger(Variant impl) : impl_(std::move(impl)) {}

				template<typename DataT>
				void log(Index iter, const std::vector<DataT>& perDOFAbs, DataT flopsThisIter = DataT(0)) const {
					std::visit([&](const auto& l){ l.log(iter, perDOFAbs, flopsThisIter); }, impl_);
				}

				template<typename DataT>
				std::vector<DataT> computePerDOFNorms(const DataT* r, Index totalSize) const {
					return std::visit([&](const auto& l){ return l.template computePerDOFNorms<DataT>(r, totalSize); }, impl_);
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
} // namespace residuum

#endif
