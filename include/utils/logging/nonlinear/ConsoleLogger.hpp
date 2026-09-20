#ifndef PDESOLVER_UTILS_LOGGING_NONLINEAR_CONSOLELOGGER_HPP
#define PDESOLVER_UTILS_LOGGING_NONLINEAR_CONSOLELOGGER_HPP

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/Types.hpp"

#include "utils/logging/core/CsvWriter.hpp"
#include "utils/logging/core/TeeStreamBuf.hpp"

namespace pdesolver {
	namespace utils {
		namespace logging {
			namespace nonlinear {

				struct ConsoleLogger {

					// config
					std::string equationName;
					std::string solverName;

					// consoleEnabled=false with textFilePath empty too means this logger prints
					// nothing at all (equivalent to NullLogger) -- callers should use NullLogger
					// directly in that case; this constructor doesn't special-case it.
					explicit ConsoleLogger(std::string equationNameIn, std::string solverNameIn, bool consoleEnabled = true, const std::string& textFilePath = "", const std::string& csvFilePath = "") :
						equationName(std::move(equationNameIn)), solverName(std::move(solverNameIn)) {

						// see solver::ConsoleLogger's identical note: teeBuf_/textFile_/out_ are
						// heap-allocated so a move of this struct (temporary -> variant) doesn't
						// dangle the pointer out_ holds into teeBuf_.
						teeBuf_ = std::make_unique<TeeStreamBuf>();

						if (consoleEnabled) teeBuf_->addTarget(std::cout.rdbuf());

						if (!textFilePath.empty()) {
							textFile_ = std::make_unique<std::ofstream>(textFilePath, std::ios::app);
							if (!textFile_->is_open()) {
								throw std::runtime_error("nonlinear::ConsoleLogger: could not open '" + textFilePath + "' for writing");
							}
							teeBuf_->addTarget(textFile_->rdbuf());
						}

						out_ = std::make_unique<std::ostream>(teeBuf_.get());

						csv_ = std::make_unique<CsvWriter>(csvFilePath, std::vector<std::string>{"iter", "residual_norm", "residual_rel", "elapsed_s"});

						startTime_ = std::chrono::steady_clock::now();

					}

					void log(Index iter, Real residualNorm, Real residualRel) const {

						const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime_).count();

						// CSV row -- every outer iteration
						if (csv_->enabled()) {
							csv_->writeRow({static_cast<Real>(iter), residualNorm, residualRel, static_cast<Real>(elapsed)});
						}

						lastIter_ = iter;
						lastResidualRel_ = residualRel;

						if (printHeader_) {
							printBanner();
							printColumnHeader();
							printHeader_ = false;
						}

						*out_ << "  " << std::left << std::setw(6) << iter << "  ";
						*out_ << std::scientific << std::setprecision(4);
						*out_ << std::setw(14) << residualNorm << "  " << std::setw(14) << residualRel << "  ";
						*out_ << std::fixed << std::setprecision(3) << elapsed << "s\n";

					}

					void summary(bool converged) const {

						const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime_).count();

						*out_ << "  [" << solverName << "] " << (converged ? "converged" : "did not converge") << ": " << lastIter_ << " iterations, res_rel " << std::scientific << std::setprecision(4) << lastResidualRel_ << ", " << std::fixed << std::setprecision(3) << elapsed << "s\n";

					}

				private:

					void printBanner() const {

						const int width = 60;

						*out_ << "\n";
						*out_ << "  " << std::string(width, '=') << "\n";
						*out_ << "  " << equationName << " - " << solverName << "\n";
						*out_ << "  " << std::string(width, '=') << "\n";

					}

					void printColumnHeader() const {

						*out_ << "  " << std::left << std::setw(6) << "Iter" << "  ";
						*out_ << std::setw(14) << "Res" << "  " << std::setw(14) << "Res[rel]" << "  ";
						*out_ << "Elapsed" << "\n";
						*out_ << "  " << std::string(60, '-') << "\n";

					}

					std::unique_ptr<TeeStreamBuf> teeBuf_;
					std::unique_ptr<std::ofstream> textFile_;
					std::unique_ptr<std::ostream> out_;
					std::unique_ptr<CsvWriter> csv_;

					mutable bool printHeader_ = true;
					mutable Index lastIter_ = 0;
					mutable Real lastResidualRel_ = Real(0);
					std::chrono::steady_clock::time_point startTime_;

				}; // struct ConsoleLogger

			} // namespace nonlinear
		} // namespace logging
	} // namespace utils
} // namespace pdesolver

#endif
