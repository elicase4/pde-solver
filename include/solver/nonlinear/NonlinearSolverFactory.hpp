#ifndef PDESOLVER_SOLVER_NONLINEAR_NONLINEARSOLVERFACTORY_HPP
#define PDESOLVER_SOLVER_NONLINEAR_NONLINEARSOLVERFACTORY_HPP

#include <memory>
#include <stdexcept>
#include <string>

#include "solver/config/NonlinearSolverConfig.hpp"
#include "solver/nonlinear/NonlinearSolverRunner.hpp"
#include "solver/stage/NonlinearCapableStage.hpp"

namespace pdesolver {
	namespace solver {
		namespace nonlinear {

			// Mirrors solver::linear::makeLinearSolverRunner's shape and throw-for-unimplemented
			// pattern. No concrete runner (NewtonRunner/PicardRunner) exists yet -- both cases
			// throw until #40 lands. Kept here now so the composition point (Stage -> Factory ->
			// Runner) is real and callable, rather than leaving HeatStage with nothing to call.
			template<stage::NonlinearCapableStage StageT, typename VectorType>
			std::unique_ptr<NonlinearSolverRunner<VectorType>> makeNonlinearSolverRunner(StageT& stage, const config::NonlinearSolverConfig& cfg, const std::string& equationLabel) {

				switch (cfg.type) {

					case config::NonlinearSolverConfig::Type::Newton:
						throw std::runtime_error("NonlinearSolverFactory[" + equationLabel + "]: Newton solver not yet implemented");

					case config::NonlinearSolverConfig::Type::Picard:
						throw std::runtime_error("NonlinearSolverFactory[" + equationLabel + "]: Picard solver not yet implemented");

				}

				throw std::runtime_error("NonlinearSolverFactory[" + equationLabel + "]: unknown nonlinear solver type");

			}

		} // namespace nonlinear
	} // namespace solver
} // namespace pdesolver

#endif
