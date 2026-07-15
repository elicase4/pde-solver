#ifndef PDESOLVER_BOUNDARYAPPLICATOR_HPP
#define PDESOLVER_BOUNDARYAPPLICATOR_HPP

#include <cstring>
#include <unordered_set>

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
				
				template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Quadrature, typename Registry>
				static void applyEssentialBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Registry& bcRegistry, const Real time, linalg::types::Vector<Real, Backend>& F);
				
				template<Index numDOFs, eval::EvalElement EvalEle, typename EvalQP, typename Quadrature, typename Registry>
				static void applyNaturalBCs(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Registry& bcRegistry, const Real time, linalg::types::Vector<Real, Backend>& F);

			}; // class BoundaryApplicator

		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#include "backend/cpu/BoundaryApplicator.tpp"
//#include "backend/cuda/BoundaryApplicator.tpp"

#endif
