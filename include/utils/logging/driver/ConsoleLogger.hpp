#ifndef PDESOLVER_UTILS_LOGGING_DRIVER_CONSOLELOGGER_HPP
#define PDESOLVER_UTILS_LOGGING_DRIVER_CONSOLELOGGER_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <utility>

#include "utils/logging/core/AnsiColor.hpp"

namespace pdesolver {
	namespace utils {
		namespace logging {
			namespace driver {

				// Lightweight narrative-event logger shared by app dispatchers (mesh, heateq)
				// and solver::driver::{Steady,Transient} -- "[tag] message", not the
				// iteration/residual reporting that utils::logging::ConsoleLogger does.
				struct ConsoleLogger {

					std::string tag;

					template<typename Args>
					void event(Args&& msg) const {
						std::cout << "[" << tag << "] " << msg << "\n";
					}

					template<typename Args>
					void warn(Args&& msg) const {
						std::ostringstream oss;
						oss << "[" << tag << "] " << msg;
						std::cout << colorize(oss.str(), Color::Yellow) << "\n";
					}

					template<typename Args>
					void error(Args&& msg) const {
						std::ostringstream oss;
						oss << "[" << tag << "] " << msg;
						std::cerr << colorize(oss.str(), Color::Red) << "\n";
					}

				}; // struct ConsoleLogger

			} // namespace driver
		} // namespace logging
	} // namespace utils
} // namespace pdesolver

#endif
