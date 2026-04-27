#ifndef PDESOLVER_CONSOLELOGGER_HPP
#define PDESOLVER_CONSOLELOGGER_HPP

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "core/Types.hpp"
#include "fem/dof/DOFOrdering.hpp"

namespace pdesolver {
	namespace utils {
		namespace logging {

			struct ConsoleLogger {
				
				// config
				std::string solverName;
				std::string label;
				std::vector<std::string> dofNames; // one name per field component
				Index interval; // print every N iterations
				bool printHeader;

				// Fields needed for per-DOF residual computation
				Index freeDOFsPerField;
				fem::dof::DOFOrdering dofOrdering;

				// single-field constructor
				explicit ConsoleLogger(std::string solverNameIn, std::string equationLabel = "Solver", Index reportInterval = 1) : solverName(std::move(SolverNameIn)), label(std::move(equationLabel)), dofNames({label}), interval(reportInterval), printHeader(true), freeDOFsPerField(0), dofOrdering(fem::dof::DOFOrdering::Interleaved) {}

				// multi-field constructor
				explicit ConsoleLogger(std::string solverNameIn, std::string equationLabel = "Solver", std::vector<std::string> componentNames, Index freeDOFsPerFieldIn, fem::dof::DOFOrdering ordering, Index reportInterval = 1) : solverName(std::move(SolverNameIn)), label(std::move(equationLabel)), dofNames(std::move(componentNames)), interval(reportInterval), printHeader(true), freeDOFsPerField(freeDOFsPerFieldIn), dofOrdering(ordering) {}

				// log()
				template<typename DataType>
				void log(Index iter, DataType absRes, DataType relRes = DataType(-1), const std::vector<DataType>& perDOF = {}) const;

				// computePerDOFNorms() helper
				template<typename DataType>
				std::vector<DataType> computePerDOFNorms(const DataType* r, Index totalSize) const;
					
				template<typename Args>
				void event(Args&& msg) const {
					std::cout << "[" << label << "]" << msg << "\n";
				}

			private:
				
				printBanner(std::string solverName) const;

				printColumnHeader(Index numPerDOF) const;

			}; // struct ConsoleLogger

		} // namespace logging
	} // namespace utils
} // namespace pdesolver

#endif
