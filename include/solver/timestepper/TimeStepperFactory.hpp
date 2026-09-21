#ifndef PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPERFACTORY_HPP
#define PDESOLVER_SOLVER_TIMESTEPPER_TIMESTEPPERFACTORY_HPP

#include <memory>
#include <stdexcept>
#include <string>

#include "solver/config/LoggingConfig.hpp"
#include "solver/config/TimeStepperConfig.hpp"
#include "solver/timestepper/BackwardEuler.hpp"
#include "solver/timestepper/StepSizePolicyFactory.hpp"
#include "solver/timestepper/TimeStepper.hpp"
#include "solver/timestepper/TimeStepperRunner.hpp"

#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/timestepper/Logger.hpp"

namespace pdesolver {
	namespace solver {
		namespace timestepper {

			inline utils::logging::timestepper::Logger makeTimestepperLogger(const config::TimeStepperConfig& cfg, const config::TimeStepperLoggerConfig& loggerCfg, const std::string& equationLabel, const std::string& timestepperName) {

				using LoggerT = utils::logging::timestepper::Logger;

				const bool consoleEnabled = (loggerCfg.type == config::LoggerConfig::Type::Console);
				const bool anyOutput = consoleEnabled || !loggerCfg.textFile.empty() || !loggerCfg.csvFile.empty();

				if (!anyOutput) return LoggerT(utils::logging::NullLogger{});

				return LoggerT(utils::logging::timestepper::ConsoleLogger(equationLabel, timestepperName, cfg.t0, cfg.tf, consoleEnabled, loggerCfg.textFile, loggerCfg.csvFile));

			}

			template<typename StageType>
			std::unique_ptr<TimeStepperRunner> makeTimeStepperRunner(StageType& stage, const config::TimeStepperConfig& cfg, const config::TimeStepperLoggerConfig& loggerCfg, const std::string& equationLabel) {

				switch (cfg.type) {

					case config::TimeStepperConfig::Type::BackwardEuler:
						static_assert(TimeStepper<BackwardEulerRunner<StageType>>, "BackwardEulerRunner<StageType> must satisfy the TimeStepper concept, including declaring TemporalOrder Order");
						return std::make_unique<BackwardEulerRunner<StageType>>(stage, cfg, makeStepSizePolicy(cfg.stepSize, equationLabel), makeTimestepperLogger(cfg, loggerCfg, equationLabel, "Backward Euler"));

					case config::TimeStepperConfig::Type::ForwardEuler:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: ForwardEuler not yet implemented");

					case config::TimeStepperConfig::Type::GeneralizedAlpha:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: GeneralizedAlpha not yet implemented");

					case config::TimeStepperConfig::Type::RK4:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: RK4 not yet implemented");

				}

				throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: unknown timestepper type");

			}

		} // namespace timestepper
	} // namespace solver
} // namespace pdesolver

#endif
