#ifndef PDESOLVER_FEM_BOUNDARY_BOUNDARYAPPLICATOR_HPP
#define PDESOLVER_FEM_BOUNDARY_BOUNDARYAPPLICATOR_HPP

#include <cstring>
#include <unordered_set>

#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/boundary/NaturalBoundaryRegistry.hpp"

#include "fem/eval/EvalQuadraturePointBoundary.hpp"
#include "fem/eval/EvalQuadraturePointVolume.hpp"

#include "mesh/Mesh.hpp"

#include "topology/TopologicalDOF.hpp"

#include "linalg/types/Vector.hpp"
#include "linalg/types/CSRMatrix.hpp"

namespace pdesolver {
	namespace fem {
		namespace boundary {

			namespace eval = pdesolver::fem::eval;

			template<typename Backend>
			class BoundaryApplicator {
			public:
				
				template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Model, typename FormRegistry, typename Quadrature>
				void applyEssentialBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const EssentialBoundaryRegistry& bcRegistry, const Real time, const Model& model, const FormRegistry& forms, const EvalEle& evalEle, const Quadrature& quadrature, linalg::types::Vector<Real, linalg::types::backend::CPU>& F);

				template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Quadrature>
				static void applyNaturalBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const NaturalBoundaryRegistry<EvalQP>& bcRegistry, const Real time, const EvalEle& evalEle, const Quadrature& quadrature, linalg::types::Vector<Real, Backend>& F);

			}; // class BoundaryApplicator

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#include "backend/cpu/BoundaryApplicator.tpp"
//#include "backend/cuda/BoundaryApplicator.tpp"

#endif
