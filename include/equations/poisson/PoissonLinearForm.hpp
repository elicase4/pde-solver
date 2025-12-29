#ifndef POISSON_BILINEARFORM_HPP
#define POISSON_BILINEARFORM_HPP

#include "fem/eval/ElementEval.hpp"
#include "fem/form/LinearForm.hpp"
#include "fem/geometry/JacobianTransform.hpp"

#include "apps/poisson/PoissonSourceTerm.hpp"

using namespace pdesolver;

template<Int Dim, Int NodesPerElement>
class PoissonLinearForm : public fem::form::LinearForm<Dim, NodesPerElement> {};

template<Int NodesPerElement>
class PoissonLinearForm<2, NodesPerElement> : public fem::form::LinearForm<2, NodesPerElement> {
public:
	PDE_HOST PDE_DEVICE static void computeElementVector(const fem::form::eval::ElementEval<2, NodesPerElement>& eleEval, Real* fe){
		
		// get physical coordinates
		Real x[NodesPerElement];
		fem::geometry::JacobianTransform<2, NodesPerElement>::computeForward(eleEval.nodeCoords, eleEval.N, x);

		// element vector assembly contribution
		for (Index a = 0; a < NodesPerElement; ++a){
			fe[a] += (PoissonSourceTerm<2>::value(0.0, x) * N[a]) * eleEval.detJ * eleEval.w;
		}
	}
};

#endif
