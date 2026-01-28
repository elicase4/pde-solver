#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/ElementEval.hpp"
#include "fem/geometry/JacobianTransform.hpp"

namespace pdesolver::fem::form {

	template<Int Dim>
	struct PoissonBilinearForm; // struct PoissonBilinearForm

	template<>
	struct PoissonBilinearForm<2> {

		template<Int NodesPerElement>
		PDE_HOST PDE_DEVICE static void computeElementMatrix<NodesPerElement>(const fem::eval::ElementEval<2, NodesPerElement>& evalCxt, Real* Ke){
			
			// physical gradients
			Real dNdx[NodesPerElement];
			Real dNdy[NodesPerElement];

			// transform parametric gradients
			Geometry::transformGradient(evalCxt.invDetJ, evalCxt.dNdxi, evalCxt.dNdeta, dNdx, dNdy);

			// element matrix assembly contribution
			for (Index a = 0; a < NodesPerElement; ++a){
				for (Index b = 0; b < NodesPerElement; ++b){
					Ke[a * NodesPerElement + b] += (dNdx[a]*dNdx[b] + dNdy[a]*dNdy[b]) * evalCxt.detJ * evalCxt.w;
				}
			}
		}
		
		template<Int NodesPerElement>
		PDE_HOST PDE_DEVICE static void computeElementOperator<NodesPerElement>(const fem::eval::ElementEval<2, NodesPerElement>& evalCxt, const Real* Ue, Real* Oe){
			
			// physical gradients
			Real dNdx[NodesPerElement];
			Real dNdy[NodesPerElement];

			// transform parametric gradients
			Geometry::transformGradient(evalCxt.invDetJ, evalCxt.dNdxi, evalCxt.dNdeta, dNdx, dNdy);

			// element matrix assembly contribution
			for (Index a = 0; a < NodesPerElement; ++a){
				for (Index b = 0; b < NodesPerElement; ++b){
					Oe[a] += (dNdx[a]*dNdx[b] + dNdy[a]*dNdy[b]) * Ue[b] * evalCxt.detJ * evalCxt.w;
				}
			}
		}

	}; // struct PoissonBilinearForm<2>

} // namespace pdesolver::fem::form

#endif
