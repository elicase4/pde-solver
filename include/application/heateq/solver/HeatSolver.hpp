#ifndef PDESOLVER_APPLICATION_HEATEQ_HEATSOLVER_HPP
#define PDESOLVER_APPLICATION_HEATEQ_HEATSOLVER_HPP

#include "application/heateq/config/HeatConfig.hpp"

#include "core/Config.hpp"
#include "core/FEM.hpp"
#include "core/IO.hpp"
#include "core/Mesh.hpp"
#include "core/LinAlg.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"
#include "core/Utils.hpp"

#include "equations/heateq/HeatEquation.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {
			namespace solver {
			
				template<typename Backend, typename HeatEquationBundle>
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
					pdesolver::mesh::Mesh mesh_;
					pdesolver::topology::TopologicalDOF<HeatEqBundle::NumDOFs> topoDOF_;

					// fem infrastructure
					pdesolver::fem::assembly::Assembler<Backend> assembler_;
					pdesolver::fem::boundary::BoundaryApplicator<Backend> bcApplicator_;
					pdesolver::fem::boundary::EssentialBoundaryRegistry essentialBCs_;
					pdesolver::fem::boundary::NaturalBoundaryRegistry<HeatEqBundle::EvalQPBdy> naturalBCs_;

					// models
					HeatEqBundle::DefaultModel defaultModel_;

				}; // class HeatSolver

			} // namespace solver
		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#endif
