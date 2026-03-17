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
				
				BoundaryApplicator(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& bcRegistry) : mesh_(mesh), topoDOF_(topoDOF), bcRegistry_(bcRegistry) {}
				
				template<eval::EvalElement EvalEle, typename EvalQP, typename Model, typename Form, typename Quadrature>
				requires eval::EvalModel<Model, EvalQP> && eval::EvalModel<Model, EvalQP>
				void apply(const Real time, linalg::types::Vector<Real, Backend>& F) const {
					applyEssentialBCs<EvalEle, EvalQP, Model, Form>(time, F);
					applyNaturalBCs<EvalEle, EvalQP, Form, Quadrature>(time, F);
				}
				
			private:

				const mesh::Mesh& mesh_;
				const topology::TopologicalDOF& topoDOF_;
				const BoundaryRegistry& bcRegistry_;

				template<eval::EvalElement EvalEle, typename EvalQP, typename Model, typename Form>
				requires eval::EvalModel<Model, EvalQP>
				void applyEssentialBCs(const Real time, linalg::types::Vector<Real, Backend>& F);
				
				template<eval::EvalElement EvalEle, typename EvalQP, typename Form, typename Quadrature>
				requires eval::EvalQuadraturePoint<EvalQP, EvalEle>
				void applyNaturalBCs(const Real time, linalg::types::Vector<Real, Backend>& F);

			}; // class BoundaryApplicator

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

//#include "backend/cpu/BoundaryApplicator.tpp"
//#include "backend/cuda/BoundaryApplicator.tpp"

#endif
