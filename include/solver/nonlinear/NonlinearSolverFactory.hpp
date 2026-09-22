#ifndef RESIDUUM_SOLVER_NONLINEAR_NONLINEARSOLVERFACTORY_HPP
#define RESIDUUM_SOLVER_NONLINEAR_NONLINEARSOLVERFACTORY_HPP

#include <memory>
#include <stdexcept>
#include <string>

#include "solver/config/LoggingConfig.hpp"
#include "solver/config/NonlinearSolverConfig.hpp"
#include "solver/nonlinear/Newton.hpp"
#include "solver/nonlinear/NonlinearSolverRunner.hpp"
#include "solver/stage/NonlinearCapableStage.hpp"

#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/nonlinear/Logger.hpp"

namespace residuum {
	namespace solver {
		namespace nonlinear {

			inline utils::logging::nonlinear::Logger makeNonlinearLogger(const config::NonlinearLoggerConfig& loggerCfg, const std::string& equationLabel, const std::string& solverName) {

				using LoggerT = utils::logging::nonlinear::Logger;

				const bool consoleEnabled = (loggerCfg.type == config::LoggerConfig::Type::Console);
				const bool anyOutput = consoleEnabled || !loggerCfg.textFile.empty() || !loggerCfg.csvFile.empty();

				if (!anyOutput) return LoggerT(utils::logging::NullLogger{});

				return LoggerT(utils::logging::nonlinear::ConsoleLogger(equationLabel, solverName, consoleEnabled, loggerCfg.textFile, loggerCfg.csvFile));

			}

			template<stage::NonlinearCapableStage StageT, typename VectorT>
			std::unique_ptr<NonlinearSolverRunner<VectorT>> makeNonlinearSolverRunner(StageT& stage, const config::NonlinearSolverConfig& cfg, const config::NonlinearLoggerConfig& loggerCfg, const std::string& equationLabel) {

				switch (cfg.type) {

					case config::NonlinearSolverConfig::Type::Newton:
						return std::make_unique<NewtonRunner<StageT, VectorT>>(stage, cfg, makeNonlinearLogger(loggerCfg, equationLabel, "Newton"));

					case config::NonlinearSolverConfig::Type::Picard:
						throw std::runtime_error("NonlinearSolverFactory[" + equationLabel + "]: Picard solver not yet implemented");

				}

				throw std::runtime_error("NonlinearSolverFactory[" + equationLabel + "]: unknown nonlinear solver type");

			}

		} // namespace nonlinear
	} // namespace solver
} // namespace residuum

#endif
