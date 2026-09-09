#ifndef PDESOLVER_CONSOLELOGGER_HPP
#define PDESOLVER_CONSOLELOGGER_HPP

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

namespace pdesolver {
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

				// consoleEnabled=false with textFilePath empty too means this logger prints
				// nothing at all (equivalent to NullLogger) -- callers should use NullLogger
				// directly in that case; this constructor doesn't special-case it.
				explicit ConsoleLogger(std::string equationNameIn, std::string solverNameIn, std::string preconditionerNameIn, std::vector<std::string> dofNamesIn, std::vector<std::pair<std::string, std::string>> extraParamsIn = {}, Index freeDOFsPerFieldIn = 0, fem::dof::DOFOrdering ordering = fem::dof::DOFOrdering::Interleaved, Index reportInterval = 1, bool consoleEnabled = true, const std::string& textFilePath = "", const std::string& csvFilePath = "") :
					equationName(std::move(equationNameIn)), solverName(std::move(solverNameIn)), preconditionerName(std::move(preconditionerNameIn)), dofNames(std::move(dofNamesIn)), extraParams(std::move(extraParamsIn)), interval(reportInterval), printHeader(true), freeDOFsPerField(freeDOFsPerFieldIn), dofOrdering(ordering) {

					// teeBuf_/textFile_/out_ are all heap-allocated (not direct value members) --
					// ConsoleLogger gets MOVED (temporary -> variant -> CGRunner), and out_ holds
					// a raw pointer into teeBuf_ (and teeBuf_'s target list holds one into
					// textFile_'s streambuf); a direct value member's ADDRESS changes on move,
					// which would silently leave that pointer dangling post-move. A unique_ptr's
					// pointee address is stable across moving the unique_ptr itself.
					teeBuf_ = std::make_unique<TeeStreamBuf>();

					if (consoleEnabled) teeBuf_->addTarget(std::cout.rdbuf());

					if (!textFilePath.empty()) {
						// append, not truncate -- see driver::ConsoleLogger's identical note;
						// a future multiphysics run could construct more than one solver logger
						// against the same configured path within one process.
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

				template<typename DataType>
				void log(Index iter, const std::vector<DataType>& perDOFAbs, DataType flopsThisIter = DataType(0)) const;

				template<typename DataType>
				std::vector<DataType> computePerDOFNorms(const DataType* r, Index totalSize) const;

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
} // namespace pdesolver

#include "ConsoleLogger.tpp"

#endif
