#ifndef PDESOLVER_UTILS_LOGGING_TIMESTEPPER_LOGGER_HPP
#define PDESOLVER_UTILS_LOGGING_TIMESTEPPER_LOGGER_HPP

#include <utility>
#include <variant>

#include "core/Types.hpp"
#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/timestepper/ConsoleLogger.hpp"

namespace pdesolver {
	namespace utils {
		namespace logging {
			namespace timestepper {

				class Logger {
				public:

					using Variant = std::variant<ConsoleLogger, NullLogger>;

					explicit Logger(Variant impl) : impl_(std::move(impl)) {}

					void log(Index step, Real time, Real dt, Index attempts, Real residualNorm) const {
						std::visit([&](const auto& l){ l.log(step, time, dt, attempts, residualNorm); }, impl_);
					}

					void summary(bool converged) const {
						std::visit([&](const auto& l){ l.summary(converged); }, impl_);
					}

				private:
					Variant impl_;

				}; // class Logger

			} // namespace timestepper
		} // namespace logging
	} // namespace utils
} // namespace pdesolver

#endif
