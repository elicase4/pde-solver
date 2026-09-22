#ifndef RESIDUUM_UTILS_LOGGING_SOLVER_CONSOLELOGGER_HPP
#define RESIDUUM_UTILS_LOGGING_SOLVER_CONSOLELOGGER_HPP

#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "core/Types.hpp"
#include "fem/dof/DOFOrdering.hpp"

#include "utils/logging/core/CsvWriter.hpp"
#include "utils/logging/core/TeeStreamBuf.hpp"

namespace residuum {
	namespace utils {
		namespace logging {

			struct ConsoleLogger {

				// config
				std::string equationName;
				std::string solverName;
				std::string preconditionerName;
				std::vector<std::string> dofNames;
				std::vector<std::pair<std::string, std::string>> extraParams;
				Index interval;
				bool printHeader;

				// Fields needed for per-DOF residual computation
				Index freeDOFsPerField;
				fem::dof::DOFOrdering dofOrdering;

				// consoleEnabled=false with an empty textFilePath prints nothing; use NullLogger directly for that case instead
				explicit ConsoleLogger(std::string equationNameIn, std::string solverNameIn, std::string preconditionerNameIn, std::vector<std::string> dofNamesIn, std::vector<std::pair<std::string, std::string>> extraParamsIn = {}, Index freeDOFsPerFieldIn = 0, fem::dof::DOFOrdering ordering = fem::dof::DOFOrdering::Interleaved, Index reportInterval = 1, bool consoleEnabled = true, const std::string& textFilePath = "", const std::string& csvFilePath = "") :
					equationName(std::move(equationNameIn)), solverName(std::move(solverNameIn)), preconditionerName(std::move(preconditionerNameIn)), dofNames(std::move(dofNamesIn)), extraParams(std::move(extraParamsIn)), interval(reportInterval), printHeader(true), freeDOFsPerField(freeDOFsPerFieldIn), dofOrdering(ordering) {

					// heap-allocated so out_'s raw pointer into teeBuf_ stays valid when ConsoleLogger is moved
					teeBuf_ = std::make_unique<TeeStreamBuf>();

					if (consoleEnabled) teeBuf_->addTarget(std::cout.rdbuf());

					if (!textFilePath.empty()) {
						// appends rather than truncates so multiple loggers can share one configured path
						textFile_ = std::make_unique<std::ofstream>(textFilePath, std::ios::app);
						if (!textFile_->is_open()) {
							throw std::runtime_error("ConsoleLogger: could not open '" + textFilePath + "' for writing");
						}
						teeBuf_->addTarget(textFile_->rdbuf());
					}

					out_ = std::make_unique<std::ostream>(teeBuf_.get());

					std::vector<std::string> csvColumns = {"outer_tick", "iter"};
					for (const auto& name : dofNames) csvColumns.push_back("abs[" + name + "]");
					for (const auto& name : dofNames) csvColumns.push_back("rel[" + name + "]");
					csvColumns.push_back("flops");
					csvColumns.push_back("elapsed_s");
					csv_ = std::make_unique<CsvWriter>(csvFilePath, std::move(csvColumns));

				}

				template<typename DataT>
				void log(Index iter, const std::vector<DataT>& perDOFAbs, DataT flopsThisIter = DataT(0)) const;

				template<typename DataT>
				std::vector<DataT> computePerDOFNorms(const DataT* r, Index totalSize) const;

				void summary(bool converged) const;

				template<typename Args>
				void event(Args&& msg) const {
					*out_ << "[" << solverName << "] " << msg << "\n";
				}

			private:

				void printBanner() const;
				void printColumnHeader() const;

				std::unique_ptr<TeeStreamBuf> teeBuf_;
				std::unique_ptr<std::ofstream> textFile_;
				std::unique_ptr<std::ostream> out_;
				std::unique_ptr<CsvWriter> csv_;

				mutable bool baselineCaptured_ = false;
				mutable std::vector<Real> initialPerDOFNorms_;
				mutable std::vector<Real> lastRelPerDOF_;
				mutable std::chrono::steady_clock::time_point startTime_;
				mutable double totalFlops_ = 0.0;
				mutable Index lastIter_ = 0;

			}; // struct ConsoleLogger

		} // namespace logging
	} // namespace utils
} // namespace residuum

#include "ConsoleLogger.tpp"

#endif
