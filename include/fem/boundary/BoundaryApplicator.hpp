#ifndef PDESOLVER_BOUNDARYAPPLICATOR_HPP
#define PDESOLVER_BOUNDARYAPPLICATOR_HPP

#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/boundary/BoundaryRegistry.hpp"

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
				
				template<eval::EvalElement EvalEle, typename EvalQP, typename Basis, typename Form, typename Model>
				requires eval::EvalModel<Model, EvalQP>
				static void applyEssentialBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& bcRegistry, const Real time, const Model& model, const Form& form, linalg::types::Vector<Real, Backend>& F);
				
				template<eval::EvalElement EvalEle, typename EvalQP, typename Basis, typename Form, typename Quadrature>
				requires eval::EvalQuadraturePoint<EvalQP, EvalEle> && eval::
				static void applyNaturalBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& bcRegistry, const Real time, const Form& form, linalg::types::Vector<Real, Backend>& F);

			}; // class BoundaryApplicator

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

//#include "backend/cpu/BoundaryApplicator.tpp"
//#include "backend/cuda/BoundaryApplicator.tpp"

#endif
