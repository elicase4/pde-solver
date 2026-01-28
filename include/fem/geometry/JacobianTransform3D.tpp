namespace pdesolver::fem::geometry{

template<Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<3, NodesPerElement>::mapToPhysical(const Real* nodeCoords, const Real* N, Real* x){
	for (Index i = 0; i < NodesPerElement; ++i){
		x[0] += N[i] * nodeCoords[3*i];
		x[1] += N[i] * nodeCoords[3*i+1];
		x[2] += N[i] * nodeCoords[3*i+2];
	}
}

template<Int NodesPerElement>
PDE_HOST PDE_DEVICE Real JacobianTransform<3, NodesPerElement>::computeMetric(const Real* nodeCoords, const Real* dNdxi, const Real* dNdeta, const Real* dNdzeta, Real* J){
	// initialize
	J[0] = 0.0; J[1] = 0.0; J[2] = 0.0;
	J[3] = 0.0; J[4] = 0.0; J[5] = 0.0;
	J[6] = 0.0; J[7] = 0.0; J[8] = 0.0;

	// compute J entries
	for (Index i = 0; i < NodesPerElement; ++i){
		J[0] += dNdxi[i]  * nodeCoords[3*i];
		J[1] += dNdeta[i] * nodeCoords[3*i];
		J[2] += dNdzeta[i] * nodeCoords[3*i];
		J[3] += dNdxi[i]  * nodeCoords[3*i+1];
		J[4] += dNdeta[i] * nodeCoords[3*i+1];
		J[5] += dNdzeta[i] * nodeCoords[3*i+1];
		J[6] += dNdxi[i]  * nodeCoords[3*i+2];
		J[7] += dNdeta[i] * nodeCoords[3*i+2];
		J[8] += dNdzeta[i] * nodeCoords[3*i+2];
	}

	// compute det(J)
	Real detJ =  (J[0] * ( (J[4] * J[8]) - (J[5] * J[7]) ))
				-(J[1] * ( (J[3] * J[8]) - (J[5] * J[6]) ))
				+(J[2] * ( (J[3] * J[7]) - (J[4] * J[6]) ));
	
	return detJ;
}

template<Int NodesPerElement>
PDE_HOST PDE_DEVICE Real JacobianTransform<3, NodesPerElement>::invertMetric(const Real* J, const Real detJ, Real* invJ){
	
	// compute determinant inverse
	Real detInvJ = 1.0 / detJ;

	// compute invJ entries
	invJ[0] = detInvJ  * ( (J[4] * J[8]) - (J[5] * J[7]) );
	invJ[1] = -detInvJ * ( (J[3] * J[8]) - (J[5] * J[6]) );
	invJ[2] = detInvJ  * ( (J[3] * J[7]) - (J[4] * J[6]) );
	invJ[3] = -detInvJ * ( (J[1] * J[8]) - (J[2] * J[7]) );
	invJ[4] = detInvJ  * ( (J[0] * J[8]) - (J[2] * J[6]) );
	invJ[5] = -detInvJ * ( (J[0] * J[7]) - (J[1] * J[6]) );
	invJ[6] = detInvJ  * ( (J[1] * J[5]) - (J[2] * J[4]) );
	invJ[7] = -detInvJ * ( (J[0] * J[5]) - (J[2] * J[3]) );
	invJ[8] = detInvJ  * ( (J[0] * J[4]) - (J[1] * J[3]) );
	
	return detInvJ;
}

template<Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<3, NodesPerElement>::transformGradient(const Real* invJ, const Real* dNdxi, const Real* dNdeta, const Real* dNdzeta, Real* dNdx, Real* dNdy, Real* dNdz){
	
	for (Index i = 0; i < NodesPerElement; ++i){
		dNdx[i] = (invJ[0] * dNdxi[i]) + (invJ[3] * dNdeta[i]) + (invJ[6] * dNdzeta[i]);
		dNdy[i] = (invJ[1] * dNdxi[i]) + (invJ[4] * dNdeta[i]) + (invJ[7] * dNdzeta[i]);
		dNdz[i] = (invJ[2] * dNdxi[i]) + (invJ[5] * dNdeta[i]) + (invJ[8] * dNdzeta[i]);
	}
	
}

template<Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<3, NodesPerElement>::computeNormal(const Real* J, const Index* tangentID, const Real nCoeff, Real* n){
	
	n[0] =        nCoeff * ((J[3+tangentID[0]]*J[6+tangentID[1]]) - (J[3+tangentID[1]]*J[6+tangentID[0]]));
	n[1] = -1.0 * nCoeff * ((J[  tangentID[0]]*J[6+tangentID[1]]) - (J[  tangentID[1]]*J[6+tangentID[0]]));
	n[2] =        nCoeff * ((J[  tangentID[0]]*J[3+tangentID[1]]) - (J[  tangentID[1]]*J[3+tangentID[0]]));

}

} // namespace pdesolver::fem::geometry
