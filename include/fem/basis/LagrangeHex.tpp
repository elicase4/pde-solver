namespace pdesolver::fem::basis{

// Implementation: eval
template<int Px, int Py, int Pz>
PDE_HOST PDE_DEVICE void LagrangeHex<Px, Py, Pz>::eval(const Real* xi, Real* N){
	Real Nx[Px + 1];
	Real Ny[Py + 1];
	Real Nz[Pz + 1];

	BasisX::eval(xi[0], Nx);
	BasisY::eval(xi[1], Ny);
	BasisZ::eval(xi[2], Nz);

	// compute tensor product
	Index a = 0;
	for (Index k = 0; k <= Pz; ++k){
		for (Index j = 0; j <= Py; ++j){
			for (Index i = 0; i <= Px; ++i){
				N[a] = Nx[i] * Ny[j] * Nz[k];
				a++;
			}
		}
	}
}

// Implementation: evalGradient
template<int Px, int Py, int Pz>
PDE_HOST PDE_DEVICE void LagrangeHex<Px, Py, Pz>::evalGradient(const Real* xi, Real* dNdxi){
	Real Nx[Px + 1], Ny[Py + 1], Nz[Pz + 1];
	Real dNx[Px + 1], dNy[Py + 1], dNz[Pz + 1];

	BasisX::eval(xi[0], Nx);
	BasisX::evalFirstDerivative(xi[0], dNx);
	BasisY::eval(xi[1], Ny);
	BasisX::evalFirstDerivative(xi[1], dNy);
	BasisZ::eval(xi[2], Nz);
	BasisX::evalFirstDerivative(xi[2], dNz);

	// compute tensor product & chain rule
	Index a = 0;
	for (Index k = 0; k <= Pz; ++k){
		for (Index j = 0; j <= Py; ++j){
			for (Index i = 0; i <= Px; ++i){
				dNdxi[a*ParametricDim    ] = dNx[i] * Ny[j] * Nz[k];
				dNdxi[a*ParametricDim + 1] = Nx[i] * dNy[j] * Nz[k];
				dNdxi[a*ParametricDim + 2] = Nx[i] * Ny[j] * dNz[k];
				a++;
			}
		}
	}
}

// Implementation: evalHessian
template<int Px, int Py, int Pz>
PDE_HOST PDE_DEVICE void LagrangeHex<Px, Py, Pz>::evalHessian(const Real* xi, Real* d2Nd2xi){
	Real Nx[Px + 1], Ny[Py + 1], Nz[Pz + 1];
	Real dNx[Px + 1], dNy[Py + 1], dNz[Pz + 1];
	Real d2Nx[Px + 1], d2Ny[Py + 1], d2Nz[Pz + 1];

	BasisX::eval(xi[0], Nx);
	BasisX::evalFirstDerivative(xi[0], dNx);
	BasisX::evalSecondDerivative(xi[0], d2Nx);
	BasisY::eval(xi[1], Ny);
	BasisX::evalFirstDerivative(xi[1], dNy);
	BasisX::evalSecondDerivative(xi[1], d2Ny);
	BasisZ::eval(xi[2], Nz);
	BasisX::evalFirstDerivative(xi[2], dNz);
	BasisX::evalSecondDerivative(xi[2], d2Nz);

	const Index NumEntries = 6;

	// compute tensor product & chain rule
	Index a = 0;
	for (Index k = 0; k <= Pz; ++k){
		for (Index j = 0; j <= Py; ++j){
			for (Index i = 0; i <= Px; ++i){
				d2Nd2xi[a*NumEntries    ] = d2Nx[i] * Ny[j] * Nz[k];
				d2Nd2xi[a*NumEntries + 1] = Nx[i] * d2Ny[j] * Nz[k];
				d2Nd2xi[a*NumEntries + 2] = Nx[i] * Ny[j] * d2Nz[k];
				d2Nd2xi[a*NumEntries + 3] = dNx[i] * dNy[j] * Nz[k];
				d2Nd2xi[a*NumEntries + 4] = dNx[i] * Ny[j] * dNz[k];
				d2Nd2xi[a*NumEntries + 5] = Nx[i] * dNy[i] * dNz[k];
				a++;
			}
		}
	}
}

// Implementation: evalLaplacian
template<int Px, int Py, int Pz>
PDE_HOST PDE_DEVICE void LagrangeHex<Px, Py, Pz>::evalLaplacian(const Real* xi, Real* lapN){
	Real Nx[Px + 1], Ny[Py + 1], Nz[Pz + 1];
	Real d2Nx[Px + 1], d2Ny[Py + 1], d2Nz[Pz + 1];

	BasisX::eval(xi[0], Nx);
	BasisX::evalSecondDerivative(xi[0], d2Nx);
	BasisY::eval(xi[1], Ny);
	BasisX::evalSecondDerivative(xi[1], d2Ny);
	BasisZ::eval(xi[2], Nz);
	BasisX::evalSecondDerivative(xi[2], d2Nz);

	// compute tensor product
	Index a = 0;
	for (Index k = 0; k <= Pz; ++k){
		for (Index j = 0; j <= Py; ++j){
			for (Index i = 0; i <= Px; ++i){
				lapN[a] = d2Nx[i] * Ny[j] * Nz[k] + Nx[i] * d2Ny[j] * Nz[k] + Nx[i] * Ny[j] * d2Nz[k];
				a++;
			}
		}
	}
}

template<int Px, int Py, int Pz>
PDE_HOST PDE_DEVICE Real LagrangeHex<Px, Py, Pz>::getFaceTopology(const Int rngID, Index* tangentID){
	switch (rngID){
		case 0:
			tangentID[0] = 1;
			tangentID[1] = 2;
			return -1.0;
		case 1:
			tangentID[0] = 1;
			tangentID[1] = 2;
			return 1.0;
		case 2:
			tangentID[0] = 0;
			tangentID[1] = 2;
			return -1.0;
		case 3:
			tangentID[0] = 0;
			tangentID[1] = 2;
			return 1.0;
		case 4:
			tangentID[0] = 0;
			tangentID[1] = 1;
			return -1.0;
		case 5:
			tangentID[0] = 0;
			tangentID[1] = 1;
			return 1.0;
		default:
			return 0.0;
	}
}

// Implementation: nodesPerFace
template<int Px, int Py, int Pz>
PDE_HOST PDE_DEVICE Index LagrangeHex<Px, Py, Pz>::nodesPerFace(const Index faceID){
	switch (faceID){
		case 0:
			return (Pz + 1)*(Py + 1);
		case 1:
			return (Pz + 1)*(Py + 1);
		case 2:
			return (Pz + 1)*(Px + 1);
		case 3:
			return (Pz + 1)*(Px + 1);
		case 4:
			return (Py + 1)*(Px + 1);
		case 5:
			return (Py + 1)*(Px + 1);
	}
}


// Implementation: getFaceNodes
template<int Px, int Py, int Pz>
PDE_HOST PDE_DEVICE void LagrangeHex<Px, Py, Pz>::getFaceNodes(const Index faceID, Index* nodeIDs){
	switch (faceID){
		case 0:
			for (int j = 0; j < (Pz + 1); ++j) {
				for (int i = 0; i < (Py + 1); ++i)
					nodeIDs[i + j*(Py + 1)] = (Index) i*(Px + 1) + j*(Px + 1)*(Py + 1);
			}
		case 1:
			for (int j = 0; j < (Pz + 1); ++j) {
				for (int i = 0; i < (Py + 1); ++i)
					nodeIDs[i + j*(Py + 1)] = (Index) i*(Px + 1) + j*(Px + 1)*(Py + 1) + Px;
			}
			break;
		case 2:
			for (int j = 0; j < (Pz + 1); ++j) {
				for (int i = 0; i < (Px + 1); ++i)
					nodeIDs[i + j*(Px + 1)] = (Index) i + j*(Px + 1)*(Py + 1);
			}
			break;
		case 3:
			for (int j = 0; j < (Pz + 1); ++j) {
				for (int i = 0; i < (Px + 1); ++i)
					nodeIDs[i + j*(Px + 1)] = (Index) i + j*(Px + 1)*(Py + 1) + (Px + 1)*Py;
			}
			break;
		case 4:
			for (int j = 0; j < (Py + 1); ++j) {
				for (int i = 0; i < (Px + 1); ++i)
					nodeIDs[i + j*(Px + 1)] = (Index) i + j*(Px + 1);
			}
			break;
		case 5:
			for (int j = 0; j < (Py + 1); ++j) {
				for (int i = 0; i < (Px + 1); ++i)
					nodeIDs[i + j*(Px + 1)] = (Index) i + j*(Px + 1) + (Px + 1)*(Py + 1)*Pz;
			}
			break;
	}
}

} // namespace pdesolver::fem::basis
