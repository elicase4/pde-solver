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
				void log(Index iter, DataType absRes, DataType relRes = DataType(-1), const std::vector<DataType>& perDOF = {}) const {
					
					// check to print based on config
					if (interval == 0) return;
					if ((iter > 0) && ((iter % interval) != 0)) return;

					// print header
					if (printHeader && iter == 0) {
						printBanner(solverName);
						printColumnHeader(perDOF.size());
						const_cast<ConsoleLogger*>(this)->printHeader = false;
					}
			
					// set print format
					std::cout << "  " << std::left <<std::setw(6) << iter << "  ";
					std::cout << std::scientific << std::setprecision(4);
					std::cout << std::setw(14) << absRes << "  ";

					// print total residual
					if (relRes >= DataType(0)) {
						std::cout << std::setw(14) << relRes << "  ";
					} else {
						std::cout << std::setw(14) << "---" << "  ";
					}

					// print residual per dof
					for (const auto& r : perDOF) {
						std::cout << std::setw(14) << r << "  ";
					}

					std::cout << "\n";

				}

				// computePerDOFNorms() helper
				template<typename DataType>
				std::vector<DataType> computePerDOFNorms(const DataType* r, Index totalSize) const {
					
					// configure norms vecor setup
					const Index D = dofNames.size();
					if ((D <= 1) || (freeDOFsPerField == 0)) return {};
					
					// norms vector
					std::vector<DataType> norms(D, DataType(0));

					// split total residual norm to each dof component
					if (dofOrdering == fem::dof::DOFOrdering::Interleaved) {
						
						// r[i] belongs to component (i % D)
						for (Index i = 0; i < totalSize; ++i) {
							Index comp = (i % D);
							norms[comp] += r[i] * r[i];
						}
					
					} else {
						
						// r[c*freeDOFsPerField .. (c+1)*freeDOFsPerField]
						for (Index c = 0; c < D; ++c){
							Index start = c * freeDOFsPerField;
							Index end = start + freeDOFsPerField;
							if (end > totalSize) end = totalSize;
							for (Index i = start; i < end; ++i){
								norms[c] += r[i] * r[i];
							}
						}

					}
					
					// take 2-norm of each component
					for (auto& v : norms) v = std::sqrt(v);

					return norms;

				}

				template<typename Args>
				void event(Args&& msg) const {
					std::cout << "[" << label << "]" << msg << "\n";
				}

			private:
				
				printBanner(std::string solverName) const {

					const int width = 60;
					
					std::cout << "\n";
					std::cout << "  " << std::string(width, '=') << "\n";
					std::cout << "  " << label << " — " << solverName << "\n";
					
					if (dofNames.size() > 1) {
						std::cout << "  Components:";
						for (const auto& n : dofNames) std::cout << "  " << n;
						std::cout << "\n";
						std::cout << "  DOF ordering: "
								  << (dofOrdering == fem::dof::DOFOrdering::Interleaved
									  ? "Interleaved (node-major)"
									  : "Block (field-major)")
								  << "\n";
					}
					
					std::cout << "  " << std::string(width, '=') << "\n";
				
				}

				printColumnHeader(Index numPerDOF) const {
					
					std::cout << "  " << std::left
							  << std::setw(6)  << "Iter" << "  "
							  << std::setw(14) << "Abs Residual" << "  "
							  << std::setw(14) << "Rel Residual" << "  ";
					
					for (Index d = 0; d < numPerDOF && d < static_cast<Index>(dofNames.size()); ++d) {
						std::cout << std::setw(14) << ("Res[" + dofNames[d] + "]") << "  ";
					}
					
					std::cout << "\n";
					std::cout << "  " << std::string(std::max<Index>(48, 18 + 16*numPerDOF), '-') << "\n";
				
				}

			}; // struct ConsoleLogger

		} // namespace logging
	} // namespace utils
} // namespace pdesolver

#endif
