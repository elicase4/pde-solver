#ifndef PDESOLVER_APPLICATION_HEATEQ_HEATSTAGE_HPP
#define PDESOLVER_APPLICATION_HEATEQ_HEATSTAGE_HPP

#include <memory>
#include <string>
#include <filesystem>
#include <stdexcept>

#include "application/heateq/HeatConfig.hpp"

#include "fem/assembly/Assembler.hpp"
#include "fem/boundary/BoundaryApplicator.hpp"
#include "fem/dof/DOFOrdering.hpp"
#include "fem/form/FormRegistry.hpp"

#include "io/MeshIO.hpp"
#include "io/FieldIO.hpp"

#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"

#include "linalg/operator/CSROperator.hpp"
#include "linalg/solver/iterative/cg/Solver.hpp"
#include "linalg/solver/iterative/cg/Config.hpp"
#include "linalg/solver/iterative/cg/Workspace.hpp"
#include "linalg/solver/base/SolverReport.hpp"
#include "linalg/solver/preconditioner/Identity.hpp"

#include "mesh/Mesh.hpp"

#include "topology/TopologicalDOF.hpp"

#include "utils/logging/solver/ConsoleLogger.hpp"
#include "utils/expression/ScalarExpression.hpp"
#include "utils/expression/VectorExpression.hpp"
#include "utils/expression/TensorExpression.hpp"

namespace pdesolver {
	namespace application {
		namespace heateq {

			template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
			class HeatStage {
			public:

				using HeatEqBundle = equations::HeatEquation<BasisType::SpatialDim, BasisType, QuadratureVolumeType, QuadratureBoundaryType>;

				using VectorT = linalg::types::Vector<Real, BackendType>;
				using MatrixT = linalg::types::CSRMatrix<Real, BackendType>;

				explicit HeatStage(const HeatConfig& config);

				void initialize();
				
				void assemble();
				
				bool solve();
				
				void finalize();

				const VectorT& solution() const { return U_; }

			private:

				HeatConfig config_;
				
				mesh::Mesh mesh_;
				std::unique_ptr<topology::TopologicalDOF<HeatEqBundle::NumDOFs>> topoDOF_;

				fem::boundary::BoundaryRegistry bcRegistry_;

				fem::assembly::Assembler<BackendType> assembler_;
				fem::assembly::BoundaryApplicator<BackendType> bcApplicator_;

				MatrixT K_;
				VectorT F_;
				VectorT U_;

				typename HeatEqBundle::ConstantConductivityModel conductivityModel_;
				typename HeatEqBundle::DefaultModel defaultModel_;

			}; // class HeatStage

		} // namespace heateq
	} // namespace application
} // namespace pdesolver

#include "application/heateq/HeatStage.tpp"

#endif
