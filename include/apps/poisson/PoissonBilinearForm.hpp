#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "fem/eval/ElementEval.hpp"
#include "fem/form/BilinearForm.hpp"
#include "fem/geometry/JacobianTransform.hpp"

#include "fem/basis/LagrangeQuad.hpp"
#include "fem/quadrature/GaussQuadratureQuad.hpp"

using namespace pdesolver;

template<Int Dim, Int NodesPerElement>
class PoissonBilinearForm : public fem::form::BilinearForm<Dim, NodesPerElement> {};

template<Int NodesPerElement>
class PoissonBilinearForm<2> : public fem::form::BilinearForm<2, NodesPerElement> {
public:
	HOST_DEVICE static void computeElementMatrix(const fem::form::eval::ElementEval<2, NodesPerElement>& eleEval, Real* Ke){
		
		// Compute nodes in a single parametric direction
		Int NodesPerElementDir = NodesPerElement / 2;

		// get weights and coordinates
		fem::quadrature::GaussQuadratureQuad<NodesPerElementDir,NodesPerElementDir>::getPoints(eleEval.xi);
		fem::quadrature::GaussQuadratureQuad<NodesPerElementDir,NodesPerElementDir>::getWeights(eleEval.w);
		
		// compute parametric gradients
		fem::basis::LagrangeQuad<1,1>::evalGradient(eleEval.xi, eleEval.dNdxi, eleEval.dNdeta);

		// Jacobian Info
		eleEval.detJ = fem::geometry::JacobianTransform<2, NodesPerElement>::computeForward(eleEval.xi, eleEval.dNdxi, eleEval.dNdeta, eleEval.J);
		eleEval.invDetJ = fem::geometry::JacobianTransform<2, NodesPerElement>::invertForward(eleEval.J, eleEval.detJ, eleEval.invJ);

		// physical gradients
		Real dNdx[NodesPerElement];
		Real dNdy[NodesPerElement];

		// transform parametric gradients
		fem::geometry::JacobianTransform<2, NodesPerElement>::transformGradient(eleEval.invDetJ, eleEval.dNdxi, eleEval.dNdeta, dNdx, dNdy);

		// TODO: Check matrix assembly loop

		// element matrix assembly loop
		for (Index a = 0; a < NodesPerElementDir; ++a){
			for (Index b = 0; b < NodesPerElementDir; ++b){
				Ke[a * NodesPerElementDir + b] += (dNdx[a]*dNdx[b] + dNdy[a]*dNdy[b]) * eleEval.detJ * w[a * NodesPerElementDir + b];
			}
		}
	}
};

#endif
