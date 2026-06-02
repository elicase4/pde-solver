#ifndef PDESOLVER_SOLVER_DRIVER_STEADY_HPP
#define PDESOLVER_SOLVER_DRIVER_STEADY_HPP

#include "core/Types.hpp"
#include "core/Fem.hpp"
#include "core/LinAlg.hpp"

namespace pdesolver {
	namespace solver {
		namespace driver {

			template<typename Backend, typename EvalElement, typename EvalQPVolume, typename EvalQPBoundary, typename MatrixForm, typename VolumeForm, typename BoundaryForm, typename VolumeModel, typename QuadratureVolume, typename QuadratureBoundary, typename BCFunction, typename SolverConfig, typename SolverWorkspace, typename SolverReport, typename SolverLogger = utils::logging::NullLogger>
			class Steady {
			public:

				struct RunResult {
					SolverReport solverReport_;
					bool success = false;
				}; // struct RunResult

				Steady(const mesh::Mesh& mesh, topology::TopologicalDOF& topoDOF, fem::boundary::BoundaryRegistry& bcRegistry, const VolumeModel& model, const MatrixForm& matrixForm, const VolumeForm& volumeForm, const BoundaryForm& boundaryForm, const SolverConfig solverConfig) : mesh_(mesh), topoDOF_(topoDOF), bcRegistry_(bcRegistry), model_(model), matrixForm_(matrixForm), volumeForm_(volumeForm), boundaryForm_(boundaryForm), solverConfig_(solverConfig) {}
				
				// run solution
				RunResult run(Real time = 0.0);

				// access solution
				const VectorType& solution() const;

			private:
				
				const mesh::Mesh& mesh_;
				topology::TopologicalDOF& topoDOF_;
				fem::boundary::BoundaryRegistry& bcRegistry_;
				
				VolumeModel model_;
				
				MatrixForm matrixForm_;
				VolumeForm volumeForm_;
				BoundaryForm boundaryForm_;

				SolverConfig solverConfig_;
				SolverWorkspace solverWorkspace_;
				SolverLogger solverLogger_;

				fem::assembly::Assembler<Backend> assembler_;
				fem::boundary::BoundaryApplicator<Backend> bcApplicator_;
				linalg::types::Vector<Real, Backend> U_, F_;
				linalg::types::CSRMatrix<Real, Backend> K_;
				SolverWorkspace solverWorkspace_;

			}; // class SteadyDriver

		} // namespace driver
	} // namespace solver
} // namespace pdesolver

#endif
