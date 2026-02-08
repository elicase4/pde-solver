namespace pdesolver::fem::basis{

// Implementation: eval
template<int Px, int Py>
PDE_HOST PDE_DEVICE void LagrangeQuad<Px, Py>::eval(const Real* xi, Real* N){
	Real Nx[Px + 1];
	Real Ny[Py + 1];
	
	BasisX::eval(xi[0], Nx);
	BasisY::eval(xi[1], Ny);

	// compute tensor product
	Index a = 0;
	for (Index j = 0; j <= Py; ++j){
		for (Index i = 0; i <= Px; ++i){
			N[a] = Nx[i] * Ny[j];
			a++;
		}
	}
}

// Implementation: evalGradient
template<int Px, int Py>
PDE_HOST PDE_DEVICE void LagrangeQuad<Px, Py>::evalGradient(const Real* xi, Real* dNdxi){
	
	Index pD = 2;
	
	Real Nx[Px + 1], Ny[Py + 1];
	Real dNx[Px + 1], dNy[Py + 1];

	BasisX::eval(xi[0], Nx);
	BasisX::evalFirstDerivative(xi[0], dNx);
	BasisY::eval(xi[1], Ny);
	BasisY::evalFirstDerivative(xi[1], dNy);

	// evaluate tensor product & chain rule
	Index a = 0;
	for (Index j = 0; j <= Py; ++j){
		for (Index i = 0; i <= Px; ++i){
			dNdxi[a*pD    ] = dNx[i] * Ny[j];
			dNdxi[a*pD + 1] = Nx[i] * dNy[j];
			a++;
		}
	}
}

// Implementation: evalHessian
template<int Px, int Py>
PDE_HOST PDE_DEVICE void LagrangeQuad<Px, Py>::evalHessian(const Real* xi, Real* d2Nd2xi){
	
	Real Nx[Px + 1], Ny[Py + 1];
	Real dNx[Px + 1], dNy[Py + 1];
	Real d2Nx[Px + 1], d2Ny[Py + 1];
	
	BasisX::eval(xi[0], Nx);
	BasisX::evalFirstDerivative(xi[0], dNx);
	BasisX::evalSecondDerivative(xi[0], d2Nx);
	BasisY::eval(xi[1], Ny);
	BasisX::evalFirstDerivative(xi[1], dNy);
	BasisY::evalSecondDerivative(xi[1], d2Ny);

	const Index NumEntries = 3;

	// evaluate tensor product & chain rule
	Index a = 0;
	for (Index j = 0; j <= Py; ++j){
		for (Index i = 0; i <= Px; ++ i){
			d2Nd2xi[a*NumEntries    ] = d2Nx[i] * Ny[j];
			d2Nd2xi[a*NumEntries + 1] = Nx[i] * d2Ny[j];
			d2Nd2xi[a*NumEntries + 2] = dNx[i] * dNy[j];
			a++;
		}
	}
}

// Implementation: evalLaplacian
template<int Px, int Py>
PDE_HOST PDE_DEVICE void LagrangeQuad<Px, Py>::evalLaplacian(const Real* xi, Real* lapN){
	
	Real Nx[Px + 1], Ny[Py + 1];
	Real d2Nx[Px + 1], d2Ny[Py + 1];
	
	BasisX::eval(xi[0], Nx);
	BasisX::evalSecondDerivative(xi[0], d2Nx);
	BasisY::eval(xi[1], Ny);
	BasisY::evalSecondDerivative(xi[1], d2Ny);

	// evaluate tensor product & chain rule
	Index a = 0;
	for (Index j = 0; j <= Py; ++j){
		for (Index i = 0; i <= Px; ++ i){
			lapN[a] = d2Nx[i] * Ny[j] + Nx[i] * d2Ny[j];
			a++;
		}
	}
}

// implementation: getFaceTopology 
template<int Px, int Py>
PDE_HOST PDE_DEVICE Real LagrangeQuad<Px, Py>::getFaceTopology(const Int faceID, Index* tangentID){
	switch (faceID){
		case 0:
			tangentID[0] = 0;
			tangentID[1] = 1;
			return -1.0;
		case 1:
			tangentID[0] = 0;
			tangentID[1] = 1;
			return 1.0;
		case 2:
			tangentID[0] = 1;
			tangentID[1] = 0;
			return -1.0;
		case 3:
			tangentID[0] = 1;
			tangentID[1] = 0;
			return 1.0;
		default:
			return 0.0;
	}
}

// implemntation: nodesPerFace
template<int Px, int Py>
PDE_HOST PDE_DEVICE Index LagrangeQuad<Px, Py>::nodesPerFace(const Index FaceID){
	switch (FaceID) {
		case 0:
			return (Py + 1);
		case 1:
			return (Py + 1);
		case 2:
			return (Px + 1);
		case 3:
			return (Px + 1);
	}
}

// implemntation: getFaceNodes
template<int Px, int Py>
PDE_HOST PDE_DEVICE void LagrangeQuad<Px, Py>::getFaceNodes(const Index FaceID, Index* nodeIDs){
	switch (FaceID) {
		case 0:
			for (int i = 0; i < (Py + 1); ++i)
				nodeIDs[i] = (Index) i*(Px + 1);
			break;
		case 1:
			for (int i = 0; i < (Py + 1); ++i)
				nodeIDs[i] = (Index) i*(Px + 1) + Px;
			break;
		case 2:
			for (int i = 0; i < (Px + 1); ++i)
				nodeIDs[i] = (Index) i;
			break;
		case 3:
			for (int i = 0; i < (Px + 1); ++i)
				nodeIDs[i] = (Index) i + (Px + 1)*Py;
			break;
	}
}

} // namespace pdesolver::fem::basis
