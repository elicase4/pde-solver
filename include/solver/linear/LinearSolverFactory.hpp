#ifndef PDESOLVER_SOLVER_LINEAR_LINEARSOLVERFACTORY_HPP
#define PDESOLVER_SOLVER_LINEAR_LINEARSOLVERFACTORY_HPP

#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/Types.hpp"

#include "solver/config/LinearSolverConfig.hpp"
#include "solver/config/LoggingConfig.hpp"

#include "linalg/solver/base/LinearSolverRunner.hpp"
#include "linalg/solver/iterative/cg/CGRunner.hpp"
#include "linalg/solver/iterative/cg/Config.hpp"

#include "linalg/solver/preconditioner/Identity.hpp"

#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/solver/ConsoleLogger.hpp"
#include "utils/logging/solver/Logger.hpp"

namespace pdesolver {
	namespace solver {
		namespace linear {

			template<typename OperatorType, typename VectorType>
			std::unique_ptr<linalg::solver::LinearSolverRunner<VectorType>> makeLinearSolverRunner(const OperatorType& op, Index n, const config::LinearSolverConfig& cfg, const config::SolverLoggerConfig& loggerCfg, const std::string& equationName, const std::vector<std::string>& dofNames) {

				if (cfg.preconditioner.type != config::PreconditionerConfig::Type::Identity) {
					throw std::runtime_error("LinearSolverFactory: only the identity preconditioner is implemented so far");
				}

				const std::string preconditionerName = "Identity";

				using PreconditionerT = linalg::solver::preconditioner::Identity<VectorType>;
				using LoggerT = utils::logging::SolverLogger;

				switch (cfg.type) {

					case config::LinearSolverConfig::Type::CG: {

						using CGConfigT = linalg::solver::iterative::cg::Config<VectorType>;
						CGConfigT cgCfg;
						cgCfg.tol = static_cast<typename CGConfigT::DataType>(cfg.tolerance);
						cgCfg.maxIters = cfg.maxIterations;

						std::ostringstream tolStream;
						tolStream << std::scientific << std::setprecision(1) << cfg.tolerance;

						std::vector<std::pair<std::string, std::string>> extraParams = {
							{"Maximum Iterations", std::to_string(cfg.maxIterations)},
							{"Tolerance", tolStream.str() + " (" + (cgCfg.tolType == linalg::solver::iterative::cg::ToleranceType::Relative ? "Relative" : "Absolute") + ")"}
						};

						const bool consoleEnabled = (loggerCfg.type == config::LoggerConfig::Type::Console);
						const bool anyOutput = consoleEnabled || !loggerCfg.textFile.empty() || !loggerCfg.csvFile.empty();

						LoggerT logger = anyOutput
							? LoggerT(utils::logging::ConsoleLogger(equationName, "PCG", preconditionerName, dofNames, extraParams, n, fem::dof::DOFOrdering::Interleaved, 1, consoleEnabled, loggerCfg.textFile, loggerCfg.csvFile))
							: LoggerT(utils::logging::NullLogger{});

						return std::make_unique<linalg::solver::iterative::cg::CGRunner<OperatorType, VectorType, PreconditionerT, LoggerT>>(op, n, cgCfg, std::move(logger));

					}

					case config::LinearSolverConfig::Type::GMRES:
						throw std::runtime_error("LinearSolverFactory: GMRES not yet implemented");

					case config::LinearSolverConfig::Type::BiCGSTAB:
						throw std::runtime_error("LinearSolverFactory: BiCGSTAB not yet implemented");

					case config::LinearSolverConfig::Type::LU:
						throw std::runtime_error("LinearSolverFactory: LU not yet implemented");

				}

				throw std::runtime_error("LinearSolverFactory: unknown linear solver type");

			}

		} // namespace linear
	} // namespace solver
} // namespace pdesolver

#endif
