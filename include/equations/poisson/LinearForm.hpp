#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "equations/poisson/PoissonSourceTerm.hpp"

#include "fem/core/Types.hpp"
#include "fem/config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"
#include "fem/geometry/JacobianTransform.hpp"

namespace pdesolver::fem::form {

	template<Int Dim>
	struct PoissonLinearForm; // struct PoissonBilinearForm
	
	template<>
	struct PoissonLinearForm<2> {

		template<Int NodesPerElement>
		PDE_HOST PDE_DEVICE static void computeElementVector<NodesPerElement>(const fem::eval::ElementEval<2, NodesPerElement>& evalCxt, Real* Fe){
			
			// get physical coordinates
			Real x[NodesPerElement];
			Geometry::mapToPhysical(evalCxt.nodeCoords, evalCxt.N, x);

			// element vector assembly contribution
			for (Index a = 0; a < NodesPerElement; ++a){
				Fe[a] += (fem::eval::PoissonSourceTerm<2>::value(0.0, x) * N[a]) * evalCxt.detJ * evalCxt.w;
			}
		}

	}; // struct PoissonLinearForm<2>

} // namespace pdesolver::fem::form
