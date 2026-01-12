#ifndef PDESOLVER_ESSENTIALBC_HPP
#define PDESOLVER_ESSENTIALBC_HPP

#include "fem/boundary/BoundaryCondition.hpp"
#include "mesh/Mesh.hpp"
#include "topology/TopologicalDOF.hpp"
#include "linalg/types/SparseMatrix.hpp"
#include "linalg/types/Vector.hpp"

namespace pdesolver {
	namespace fem {
		namespace boundary {

			template<EssentialBC BC, typename Backend>
			void evaluateLift(const BC& bc, const mesh::Mesh& mesh, const dof::TopologicalDOF& topoDOF, Int boundaryTag);
			
		} // namespace boundary
	} // namespace fem
} // namespace pdesolver

#endif
