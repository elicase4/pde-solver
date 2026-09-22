#ifndef RESIDUUM_APPLICATION_HEATEQ_SOLVER_HEATSOLVER_HPP
#define RESIDUUM_APPLICATION_HEATEQ_SOLVER_HEATSOLVER_HPP

#include "application/heateq/config/HeatConfig.hpp"

#include "core/Config.hpp"
#include "core/FEM.hpp"
#include "core/IO.hpp"
#include "core/Mesh.hpp"
#include "core/LinAlg.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"
#include "core/Utils.hpp"

#include "equation/heateq/HeatEquation.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace solver {
			
				template<typename BackendT, typename HeatEqBundleT>
				class HeatSolver {
				public:

					explicit HeatSolver(const config::HeatConfig& config);

					void run();

				private:

					// methods
					void setupMesh();
					void setupTopology();
					void setupModel();
					void setupForms();
					void setupBoundaryConditions();
					void solve();
					void writeOutput() const;
					void writeLog() const;
					
					// configuration
					config::HeatConfig config_;
					
					// discretization
					residuum::mesh::Mesh mesh_;
					residuum::topology::TopologicalDOF<HeatEqBundle::NumDOFs> topoDOF_;

					// fem infrastructure
					residuum::fem::assembly::Assembler<BackendT> assembler_;
					residuum::fem::boundary::BoundaryApplicator<BackendT> bcApplicator_;
					residuum::fem::boundary::EssentialBoundaryRegistry essentialBCs_;
					residuum::fem::boundary::NaturalBoundaryRegistry<HeatEqBundle::EvalQPBdy> naturalBCs_;

					// models
					HeatEqBundle::DefaultModel defaultModel_;

				}; // class HeatSolver

			} // namespace solver
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
