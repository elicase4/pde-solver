namespace pdesolver::fem::geometry {

template<Int SpatialDim, Int ParametricDim, Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDim, ParametricDim, NodesPerElement>::mapToPhysical(const Real* nodeCoords, const Real* N, Real* x){
	
	// initialize
	for (Index i = 0; i < SpatialDim; ++i){
		x[i] = 0.0;
	}

	// compute x entries
	for (Index a = 0; a < NodesPerElement; ++a){
		for (Index i = 0; i < SpatialDim; ++i){
			x[i] += N[a]  * nodeCoords[a*SpatialDim + i];
		}
	}

}

template<Int SpatialDim, Int ParametricDim, Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDim, ParametricDim, NodesPerElement>::computeJacobian(const Real* nodeCoords, const Real* dNdxi, Real* J){
	
	// initialize
	for (Index i = 0; i < SpatialDim*ParametricDim; ++i){
		J[i] = 0.0;
	}

	// compute J entries
	for (Index i = 0; i < SpatialDim; ++i){
		for (Index alpha = 0; alpha < ParametricDim; ++alpha){
			for (Index a = 0; a < NodesPerElement; ++a){
				J[i*ParametricDim + alpha] += dNdxi[a*ParametricDim + alpha]  * nodeCoords[a*SpatialDim + i];
			}
		}
	}

}

template<Int SpatialDim, Int ParametricDim, Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDim, ParametricDim, NodesPerElement>::computeMetric(const Real* J, Real* g){

	for (Index beta = 0; beta < ParametricDim; ++beta){
		for (Index alpha = 0; alpha < ParametricDim; ++alpha){
			
			Real sum = 0.0;
			for (Index i = 0; i < SpatialDim; ++i){
				sum += J[i*ParametricDim + alpha] * J[i*ParametricDim + beta];
			}

			g[beta*ParametricDim + alpha] = sum;
			g[alpha*ParametricDim + beta] = sum;
		}
	}

}

template<Int SpatialDim, Int ParametricDim, Int NodesPerElement>
PDE_HOST PDE_DEVICE Real JacobianTransform<SpatialDim, ParametricDim, NodesPerElement>::computeMeasure(const Real* g){

	Real detg = computeMatrixDeterminant(g);
	Real measure = sqrt(detg);

	return measure;

}

template<Int SpatialDim, Int ParametricDim, Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDim, ParametricDim, NodesPerElement>::transformGradient(const Real* J, const Real* g, const Real* dNdxi, Real* dNdx){

	// initialize
	for (Index i = 0; i < NodesPerElement*SpatialDim; ++i){
		dNdx[i] = 0.0;
	}
	
	if constexpr (ParametricDim == SpatialDim){

		Real invJ[ParametricDim*ParametricDim];
		Real detJ = computeMatrixDeterminant(J);
		invertMatrix(J, detJ, invJ);
	
		for (Index a = 0; a < NodesPerElement; ++a){
			for (Index i = 0; i < SpatialDim; ++i){
				for (Index alpha = 0; alpha < ParametricDim; ++alpha){
					dNdx[a*SpatialDim + i] += (invJ[i*ParametricDim + alpha] * dNdxi[a*ParametricDim + alpha]);
				}
			}
		}

	} else if constexpr (ParametricDim != SpatialDim){

		Real invg[ParametricDim*ParametricDim];
		Real detg = computeMatrixDeterminant(g);
		invertMatrix(g, detg, invg);
	
		for (Index a = 0; a < NodesPerElement; ++a){
			for (Index i = 0; i < SpatialDim; ++i){
				for (Index beta = 0; beta < ParametricDim; ++beta){
					for (Index alpha = 0; alpha < ParametricDim; ++alpha){
						dNdx[a*SpatialDim + i] += (J[i*ParametricDim + beta] * invg[beta*ParametricDim + alpha] * dNdxi[a*ParametricDim + alpha]);
					}
				}
			}
		}

	}
	
}

template<Int SpatialDim, Int ParametricDim, Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDim, ParametricDim, NodesPerElement>::computeNormal(const Real* J, const Index* tangentID, const Real nCoeff, Real* n){
	
	if constexpr (ParametricDim == 2 && SpatialDim == 2) {

		Real t0 = J[tangentID[0]];
		Real t1 = J[ParametricDim+tangentID[0]];
		
		n[0] = nCoeff * t1;
		n[1] = -1.0 * nCoeff * t0;
		
	} else if constexpr (ParametricDim == 2 && SpatialDim == 3){
		
		Real t0 = J[tangentID[0]];
		Real t1 = J[ParametricDim+tangentID[0]];
		Real t2 = J[2*ParametricDim+tangentID[0]];
		
		Real a0 = J[tangentID[1]];
		Real a1 = J[ParametricDim+tangentID[1]];
		Real a2 = J[2*ParametricDim+tangentID[1]];

		n[0] = nCoeff * ( ((t0*a0 + t1*a1 + t2*a2) * t0) - ((t0*t0 + t1*t1 + t2*t2) * a0) );
		n[1] = nCoeff * ( ((t0*a0 + t1*a1 + t2*a2) * t1) - ((t0*t0 + t1*t1 + t2*t2) * a1) );
		n[2] = nCoeff * ( ((t0*a0 + t1*a1 + t2*a2) * t2) - ((t0*t0 + t1*t1 + t2*t2) * a2) );

	} else if constexpr (ParametricDim == 3 && SpatialDim == 3){
		
		Real t0 = J[tangentID[0]];
		Real t1 = J[ParametricDim+tangentID[0]];
		Real t2 = J[2*ParametricDim+tangentID[0]];
		
		Real a0 = J[tangentID[1]];
		Real a1 = J[ParametricDim+tangentID[1]];
		Real a2 = J[2*ParametricDim+tangentID[1]];
		
		n[0] = nCoeff * ( (t1*a2) - (a1*t2) );
		n[1] = -1.0 * nCoeff * ( (t0*a2) - (a0*t2) );
		n[2] = nCoeff * ( (t0*a1) - (a0*t1) );

	}

}

template<Int SpatialDim, Int ParametricDim, Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDim, ParametricDim, NodesPerElement>::invertMatrix(const Real* A, const Real detA, Real* invA){
	
	Real detInvA = 1.0 / detA;

	if constexpr (ParametricDim == 2) {
		
		invA[0] = detInvA  * A[3];
		invA[1] = -detInvA * A[1];
		invA[2] = -detInvA * A[2];
		invA[3] = detInvA  * A[0];

	} else if constexpr (ParametricDim == 3){

		invA[0] = detInvA  * ( (A[4] * A[8]) - (A[5] * A[7]) );
		invA[1] = -detInvA * ( (A[3] * A[8]) - (A[5] * A[6]) );
		invA[2] = detInvA  * ( (A[3] * A[7]) - (A[4] * A[6]) );
		invA[3] = -detInvA * ( (A[1] * A[8]) - (A[2] * A[7]) );
		invA[4] = detInvA  * ( (A[0] * A[8]) - (A[2] * A[6]) );
		invA[5] = -detInvA * ( (A[0] * A[7]) - (A[1] * A[6]) );
		invA[6] = detInvA  * ( (A[1] * A[5]) - (A[2] * A[4]) );
		invA[7] = -detInvA * ( (A[0] * A[5]) - (A[2] * A[3]) );
		invA[8] = detInvA  * ( (A[0] * A[4]) - (A[1] * A[3]) );

	}

}

template<Int SpatialDim, Int ParametricDim, Int NodesPerElement>
PDE_HOST PDE_DEVICE Real JacobianTransform<SpatialDim, ParametricDim, NodesPerElement>::computeMatrixDeterminant(const Real* A){

	Real detA;

	if constexpr (ParametricDim == 2) {
		
		detA = (A[0] * A[3]) - (A[1] * A[2]);

	} else if constexpr (ParametricDim == 3){
		
		detA = (A[0] * ( (A[4] * A[8]) - (A[5] * A[7]) ))
			  -(A[1] * ( (A[3] * A[8]) - (A[5] * A[6]) ))
			  +(A[2] * ( (A[3] * A[7]) - (A[4] * A[6]) ));

	}

	return detA;

}

} // namespace pdesolver::fem::geometry
