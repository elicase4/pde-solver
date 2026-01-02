#include "JacobianTransform.hpp"

using namespace fem::geometry;

template<Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<2, NodesPerElement>::computeForward(const Real* nodeCoords, const Real* N, Real* x){
	for (Index i = 0; i < NodesPerElement; ++i){
		x[0] += N[i] * nodeCoords[2*i];
		x[1] += N[i] * nodeCoords[2*i+1];
	}
}

template<Int NodesPerElement>
PDE_HOST PDE_DEVICE Real JacobianTransform<2, NodesPerElement>::computeJacobian(const Real* nodeCoords, const Real* dNdxi, const Real* dNdeta, Real* J){
	// initialize
	J[0] = 0.0; J[1] = 0.0;
	J[2] = 0.0; J[3] = 0.0;

	// compute J entries
	for (Index i = 0; i < NodesPerElement; ++i){
		J[0] += dNdxi[i]  * nodeCoords[2*i];
		J[1] += dNdeta[i] * nodeCoords[2*i];
		J[2] += dNdxi[i]  * nodeCoords[2*i+1];
		J[3] += dNdeta[i] * nodeCoords[2*i+1];
	}

	// compute det(J)
	Real detJ = (J[0] * J[3]) - (J[1] * J[2]);
	
	return detJ;
}

template<Int NodesPerElement>
PDE_HOST PDE_DEVICE Real JacobianTransform<2, NodesPerElement>::invertJacobian(const Real* J, const Real detJ, Real* invJ){
	
	// compute determinant inverse
	Real detInvJ = 1.0 / detJ;

	// compute invJ entries
	invJ[0] = detInvJ  * J[3];
	invJ[1] = -detInvJ * J[1];
	invJ[2] = -detInvJ * J[2];
	invJ[3] = detInvJ  * J[0];
	
	return detInvJ;
}

template<Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<2, NodesPerElement>::transformGradient(const Real* invJ, const Real* dNdxi, const Real* dNdeta, Real* dNdx, Real* dNdy){
	
	for (Index i = 0; i < NodesPerElement; ++i){
		dNdx[i] = (invJ[0] * dNdxi[i]) + (invJ[2] * dNdeta[i]);
		dNdy[i] = (invJ[1] * dNdxi[i]) + (invJ[3] * dNdeta[i]);
	}
	
}

template<Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<2, NodesPerElement>::computeNormal(const Real* J, const Index* tangentID, const Real nCoeff, Real* n){
	
	n[0] =        nCoeff * J[2+tangentID[0]];
	n[1] = -1.0 * nCoeff * J[  tangentID[1]];

}
