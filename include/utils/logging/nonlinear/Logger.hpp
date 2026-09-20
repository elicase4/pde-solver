#ifndef PDESOLVER_UTILS_LOGGING_NONLINEAR_LOGGER_HPP
#define PDESOLVER_UTILS_LOGGING_NONLINEAR_LOGGER_HPP

#include <utility>
#include <variant>

#include "core/Types.hpp"
#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/nonlinear/ConsoleLogger.hpp"

namespace pdesolver {
	namespace utils {
		namespace logging {
			namespace nonlinear {

				class Logger {
				public:

					using Variant = std::variant<ConsoleLogger, NullLogger>;

					explicit Logger(Variant impl) : impl_(std::move(impl)) {}

					void log(Index iter, Real residualNorm, Real residualRel) const {
						std::visit([&](const auto& l){ l.log(iter, residualNorm, residualRel); }, impl_);
					}

					void summary(bool converged) const {
						std::visit([&](const auto& l){ l.summary(converged); }, impl_);
					}

				private:
					Variant impl_;

				}; // class Logger

			} // namespace nonlinear
		} // namespace logging
	} // namespace utils
} // namespace pdesolver

#endif
