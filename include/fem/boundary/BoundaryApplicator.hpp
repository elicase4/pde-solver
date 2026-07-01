#ifndef PDESOLVER_BOUNDARYAPPLICATOR_HPP
#define PDESOLVER_BOUNDARYAPPLICATOR_HPP

#include <cstring>

#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/boundary/BoundaryRegistry.hpp"

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
				
				template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename FormRegistry, typename Model, typename Quadrature>
				static void applyEssentialBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const BoundaryRegistry& bcRegistry, const Real time, const Model& model, const FormRegistry& forms, linalg::types::Vector<Real, Backend>& F);
				
				template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename FormRegistry, typename Quadrature>
				static void applyNaturalBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const BoundaryRegistry& bcRegistry, const Real time, const FormRegistry& forms, linalg::types::Vector<Real, Backend>& F);

			}; // class BoundaryApplicator

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#include "backend/cpu/BoundaryApplicator.tpp"
//#include "backend/cuda/BoundaryApplicator.tpp"

#endif
