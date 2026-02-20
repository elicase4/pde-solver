namespace pdesolver::fem::geometry {

template<Int SpatialDimension, Int ParametricDimension, Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension, NodesPerElement>::mapToPhysical(const Real* nodeCoords, const Real* N, Real* x){
	
	// initialize
	for (Index i = 0; i < SpatialDimension; ++i){
		x[i] = 0.0;
	}

	// compute x entries
	for (Index a = 0; a < NodesPerElement; ++a){
		for (Index i = 0; i < SpatialDimension; ++i){
			x[i] += N[a]  * nodeCoords[a*SpatialDimension + i];
		}
	}

}

template<Int SpatialDimension, Int ParametricDimension, Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension, NodesPerElement>::computeJacobian(const Real* nodeCoords, const Real* dNdxi, Real* J){
	
	// initialize
	for (Index i = 0; i < SpatialDimension*ParametricDimension; ++i){
		J[i] = 0.0;
	}

	// compute J entries
	for (Index i = 0; i < SpatialDimension; ++i){
		for (Index alpha = 0; alpha < ParametricDimension; ++alpha){
			for (Index a = 0; a < NodesPerElement; ++a){
				J[i*ParametricDimension + alpha] += dNdxi[a*ParametricDimension + alpha]  * nodeCoords[a*SpatialDimension + i];
			}
		}
	}

}

template<Int SpatialDimension, Int ParametricDimension, Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension, NodesPerElement>::computeMetric(const Real* J, Real* g){

	for (Index beta = 0; beta < ParametricDimension; ++beta){
		for (Index alpha = 0; alpha < ParametricDimension; ++alpha){
			
			Real sum = 0.0;
			for (Index i = 0; i < SpatialDimension; ++i){
				sum += J[i*ParametricDimension + alpha] * J[i*ParametricDimension + beta];
			}

			g[beta*ParametricDimension + alpha] = sum;
			g[alpha*ParametricDimension + beta] = sum;
		}
	}

}

template<Int SpatialDimension, Int ParametricDimension, Int NodesPerElement>
PDE_HOST PDE_DEVICE Real JacobianTransform<SpatialDimension, ParametricDimension, NodesPerElement>::computeMeasure(const Real* g){

	Real detg = computeMatrixDeterminant(g);
	Real measure = sqrt(detg);

	return measure;

}

template<Int SpatialDimension, Int ParametricDimension, Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension, NodesPerElement>::transformGradient(const Real* J, const Real* g, const Real* dNdxi, Real* dNdx){

	// initialize
	for (Index i = 0; i < NodesPerElement*SpatialDimension; ++i){
		dNdx[i] = 0.0;
	}
	
	if constexpr (ParametricDimension == SpatialDimension){

		Real invJ[ParametricDimension*ParametricDimension];
		Real detJ = computeMatrixDeterminant(J);
		invertMatrix(J, detJ, invJ);
	
		for (Index a = 0; a < NodesPerElement; ++a){
			for (Index i = 0; i < SpatialDimension; ++i){
				for (Index alpha = 0; alpha < ParametricDimension; ++alpha){
					dNdx[a*SpatialDimension + i] += (invJ[i*ParametricDimension + alpha] * dNdxi[a*ParametricDimension + alpha]);
				}
			}
		}

	} else if constexpr (ParametricDimension != SpatialDimension){

		Real invg[ParametricDimension*ParametricDimension];
		Real detg = computeMatrixDeterminant(g);
		invertMatrix(g, detg, invg);
	
		for (Index a = 0; a < NodesPerElement; ++a){
			for (Index i = 0; i < SpatialDimension; ++i){
				for (Index beta = 0; beta < ParametricDimension; ++beta){
					for (Index alpha = 0; alpha < ParametricDimension; ++alpha){
						dNdx[a*SpatialDimension + i] += (J[i*ParametricDimension + beta] * invg[beta*ParametricDimension + alpha] * dNdxi[a*ParametricDimension + alpha]);
					}
				}
			}
		}

	}
	
}

template<Int SpatialDimension, Int ParametricDimension, Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension, NodesPerElement>::computeNormal(const Real* J, const Index* tangentID, const Real nCoeff, Real* n){
	
	if constexpr (ParametricDimension == 2 && SpatialDimension == 2) {

		Real t0 = J[tangentID[0]];
		Real t1 = J[ParametricDimension+tangentID[0]];
		
		n[0] = nCoeff * t1;
		n[1] = -1.0 * nCoeff * t0;
		
	} else if constexpr (ParametricDimension == 2 && SpatialDimension == 3){
		
		Real t0 = J[tangentID[0]];
		Real t1 = J[ParametricDimension+tangentID[0]];
		Real t2 = J[2*ParametricDimension+tangentID[0]];
		
		Real a0 = J[tangentID[1]];
		Real a1 = J[ParametricDimension+tangentID[1]];
		Real a2 = J[2*ParametricDimension+tangentID[1]];

		n[0] = nCoeff * ( ((t0*a0 + t1*a1 + t2*a2) * t0) - ((t0*t0 + t1*t1 + t2*t2) * a0) );
		n[1] = nCoeff * ( ((t0*a0 + t1*a1 + t2*a2) * t1) - ((t0*t0 + t1*t1 + t2*t2) * a1) );
		n[2] = nCoeff * ( ((t0*a0 + t1*a1 + t2*a2) * t2) - ((t0*t0 + t1*t1 + t2*t2) * a2) );

	} else if constexpr (ParametricDimension == 3 && SpatialDimension == 3){
		
		Real t0 = J[tangentID[0]];
		Real t1 = J[ParametricDimension+tangentID[0]];
		Real t2 = J[2*ParametricDimension+tangentID[0]];
		
		Real a0 = J[tangentID[1]];
		Real a1 = J[ParametricDimension+tangentID[1]];
		Real a2 = J[2*ParametricDimension+tangentID[1]];
		
		n[0] = nCoeff * ( (t1*a2) - (a1*t2) );
		n[1] = -1.0 * nCoeff * ( (t0*a2) - (a0*t2) );
		n[2] = nCoeff * ( (t0*a1) - (a0*t1) );

	}

}

template<Int SpatialDimension, Int ParametricDimension, Int NodesPerElement>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension, NodesPerElement>::invertMatrix(const Real* A, const Real detA, Real* invA){
	
	Real detInvA = 1.0 / detA;

	if constexpr (ParametricDimension == 2) {
		
		invA[0] = detInvA  * A[3];
		invA[1] = -detInvA * A[1];
		invA[2] = -detInvA * A[2];
		invA[3] = detInvA  * A[0];

	} else if constexpr (ParametricDimension == 3){

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

template<Int SpatialDimension, Int ParametricDimension, Int NodesPerElement>
PDE_HOST PDE_DEVICE Real JacobianTransform<SpatialDimension, ParametricDimension, NodesPerElement>::computeMatrixDeterminant(const Real* A){

	Real detA;

	if constexpr (ParametricDimension == 2) {
		
		detA = (A[0] * A[3]) - (A[1] * A[2]);

	} else if constexpr (ParametricDimension == 3){
		
		detA = (A[0] * ( (A[4] * A[8]) - (A[5] * A[7]) ))
			  -(A[1] * ( (A[3] * A[8]) - (A[5] * A[6]) ))
			  +(A[2] * ( (A[3] * A[7]) - (A[4] * A[6]) ));

	}

	return detA;

}

} // namespace pdesolver::fem::geometry
