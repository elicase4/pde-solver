#ifndef PDESOLVER_UTILS_LOGGING_DRIVER_CONSOLELOGGER_HPP
#define PDESOLVER_UTILS_LOGGING_DRIVER_CONSOLELOGGER_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "utils/logging/core/AnsiColor.hpp"

namespace pdesolver {
	namespace utils {
		namespace logging {
			namespace driver {

				struct ConsoleLogger {

					std::string tag;
					bool consoleEnabled = true;

					ConsoleLogger(std::string tagIn, bool consoleEnabledIn = true, const std::string& textFilePath = "") : tag(std::move(tagIn)), consoleEnabled(consoleEnabledIn) {

						if (!textFilePath.empty()) {
							// append, not truncate -- HeatDispatcher and HeatProblem each
							// construct their own driver logger against the SAME configured
							// path; two independent ofstreams truncating the same file would
							// each wipe what the other already wrote
							textFile_.open(textFilePath, std::ios::app);
							if (!textFile_.is_open()) {
								throw std::runtime_error("driver::ConsoleLogger: could not open '" + textFilePath + "' for writing");
							}
						}

					}

					template<typename Args>
					void event(Args&& msg) const {
						std::ostringstream oss;
						oss << "[" << tag << "] " << msg;
						if (consoleEnabled) std::cout << oss.str() << "\n";
						if (textFile_.is_open()) textFile_ << oss.str() << "\n" << std::flush;
					}

					template<typename Args>
					void warn(Args&& msg) const {
						std::ostringstream oss;
						oss << "[" << tag << "] " << msg;
						if (consoleEnabled) std::cout << colorize(oss.str(), Color::Yellow) << "\n";
						if (textFile_.is_open()) textFile_ << oss.str() << "\n" << std::flush;
					}

					template<typename Args>
					void error(Args&& msg) const {
						std::ostringstream oss;
						oss << "[" << tag << "] " << msg;
						if (consoleEnabled) std::cerr << colorize(oss.str(), Color::Red) << "\n";
						if (textFile_.is_open()) textFile_ << oss.str() << "\n" << std::flush;
					}

				private:

					mutable std::ofstream textFile_;

				}; // struct ConsoleLogger

			} // namespace driver
		} // namespace logging
	} // namespace utils
} // namespace pdesolver

#endif
