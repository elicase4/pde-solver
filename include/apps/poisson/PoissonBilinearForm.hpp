#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "fem/eval/ElementEval.hpp"
#include "fem/form/BilinearForm.hpp"
#include "fem/geometry/JacobianTransform.hpp"

using namespace pdesolver;

template<Int Dim, Int NodesPerElement>
class PoissonBilinearForm : public fem::form::BilinearForm<Dim, NodesPerElement> {};

template<Int NodesPerElement>
class PoissonBilinearForm<2, NodesPerElement> : public fem::form::BilinearForm<2, NodesPerElement> {
public:
	HOST_DEVICE static void computeElementMatrix(const fem::form::eval::ElementEval<2, NodesPerElement>& eleEval, Real* Ke){
		
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
};

#endif
