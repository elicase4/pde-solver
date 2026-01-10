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
		PDE_HOST PDE_DEVICE void computeElementMatrix<NodesPerElement>(const fem::eval::ElementEval<2, NodesPerElement>& eleEval, Real* Ke){
			
			// physical gradients
			Real dNdx[NodesPerElement];
			Real dNdy[NodesPerElement];

			// transform parametric gradients
			fem::geometry::JacobianTransform<2, NodesPerElement>::transformGradient(eleEval.invDetJ, eleEval.dNdxi, eleEval.dNdeta, dNdx, dNdy);

			// element matrix assembly contribution
			for (Index a = 0; a < NodesPerElement; ++a){
				for (Index b = 0; b < NodesPerElement; ++b){
					Ke[a * NodesPerElement + b] += (dNdx[a]*dNdx[b] + dNdy[a]*dNdy[b]) * eleEval.detJ * eleEval.w;
				}
			}
		}
		
		template<Int NodesPerElement>
		PDE_HOST PDE_DEVICE void computeElementOperator<NodesPerElement>(const fem::eval::ElementEval<2, NodesPerElement>& eleEval, const Real* Ue, Real* Oe){
			
			// physical gradients
			Real dNdx[NodesPerElement];
			Real dNdy[NodesPerElement];

			// transform parametric gradients
			fem::geometry::JacobianTransform<2, NodesPerElement>::transformGradient(eleEval.invDetJ, eleEval.dNdxi, eleEval.dNdeta, dNdx, dNdy);

			// element matrix assembly contribution
			for (Index a = 0; a < NodesPerElement; ++a){
				for (Index b = 0; b < NodesPerElement; ++b){
					Oe[a] += (dNdx[a]*dNdx[b] + dNdy[a]*dNdy[b]) * Ue[b] * eleEval.detJ * eleEval.w;
				}
			}
		}

	}; // struct PoissonBilinearForm<2>

} // namespace pdesolver::fem::form

#endif
