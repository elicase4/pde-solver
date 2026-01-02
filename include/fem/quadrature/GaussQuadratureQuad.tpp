namespace pdesolver::fem::quadrature {
		
// implementation points
template<Int NumPointsX, Int NumPointsY>
PDE_HOST PDE_DEVICE void GaussQuadratureQuad<NumPointsX, NumPointsY>::getPoints(Real* xi){
	
	Real xi_xi[NumPointsX];
	Real xi_eta[NumPointsY];

	QuadX::getPoints(xi_xi);
	QuadY::getPoints(xi_eta);

	// compute tensor product
	Index q;
	for (Index j = 0; j < NumPointsY; ++j) {
		for (Index i = 0; i < NumPointsX; ++i) {
			xi[2*q]   = xi_xi[i];
			xi[2*q+1] = xi_eta[j];
			q++;
		}
	}
}

// implementation weights
template<Int NumPointsX, Int NumPointsY>
PDE_HOST PDE_DEVICE void GaussQuadratureQuad<NumPointsX, NumPointsY>::getWeights(Real* w){
	
	Real w_xi[NumPointsX];
	Real w_eta[NumPointsY];

	QuadX::getPoints(w_xi);
	QuadY::getPoints(w_eta);

	// compute tensor product
	Index q;
	for (Index j = 0; j < NumPointsY; ++j) {
		for (Index i = 0; i < NumPointsX; ++i) {
			w[q] = w_xi[i] * w_eta[j];
			q++;
		}
	}
}

} // namespace pdesolver::fem::quadrature
