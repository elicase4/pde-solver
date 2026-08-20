#ifndef PDESOLVER_CONSOLELOGGER_HPP
#define PDESOLVER_CONSOLELOGGER_HPP

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "core/Types.hpp"
#include "fem/dof/DOFOrdering.hpp"

namespace pdesolver {
	namespace utils {
		namespace logging {

			struct ConsoleLogger {

				// config
				std::string equationName;    // banner headline, e.g. "Heat Equation"
				std::string solverName;      // banner headline + [tag] bracket, e.g. "PCG"
				std::string preconditionerName;
				std::vector<std::string> dofNames; // one name per field component, e.g. {"T"}
				std::vector<std::pair<std::string, std::string>> extraParams; // extra banner lines
				Index interval; // print every N iterations
				bool printHeader;

				// Fields needed for per-DOF residual computation
				Index freeDOFsPerField;
				fem::dof::DOFOrdering dofOrdering;

				explicit ConsoleLogger(std::string equationNameIn, std::string solverNameIn, std::string preconditionerNameIn, std::vector<std::string> dofNamesIn, std::vector<std::pair<std::string, std::string>> extraParamsIn = {}, Index freeDOFsPerFieldIn = 0, fem::dof::DOFOrdering ordering = fem::dof::DOFOrdering::Interleaved, Index reportInterval = 1) : equationName(std::move(equationNameIn)), solverName(std::move(solverNameIn)), preconditionerName(std::move(preconditionerNameIn)), dofNames(std::move(dofNamesIn)), extraParams(std::move(extraParamsIn)), interval(reportInterval), printHeader(true), freeDOFsPerField(freeDOFsPerFieldIn), dofOrdering(ordering) {}

				// log() -- perDOFAbs are absolute per-DOF residual norms; the logger tracks the
				// iteration-0 baseline internally and prints relative values. flopsThisIter is
				// the (usually precomputed, constant) flop cost attributable to this iteration.
				template<typename DataType>
				void log(Index iter, const std::vector<DataType>& perDOFAbs, DataType flopsThisIter = DataType(0)) const;

				// computePerDOFNorms() helper -- absolute, not relative
				template<typename DataType>
				std::vector<DataType> computePerDOFNorms(const DataType* r, Index totalSize) const;

				// one-line closing statement -- call once after the solve loop ends
				void summary(bool converged) const;

				template<typename Args>
				void event(Args&& msg) const {
					std::cout << "[" << solverName << "] " << msg << "\n";
				}

			private:

				void printBanner() const;
				void printColumnHeader() const;

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
