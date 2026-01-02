namespace pdesolver::fem::basis{

// Implementation: eval
template <int Order>
PDE_HOST PDE_DEVICE void Lagrange1D<Order>::eval(Real xi, Real* N){
	if constexpr (Order == 1){
		N[0] = 0.5 * (1.0 - xi);
		N[1] = 0.5 * (1.0 + xi);
	}
	else if constexpr (Order == 2){
		N[0] = 0.5 * xi * (xi - 1.0);
		N[1] = (1.0 - xi) * (1.0 + xi);
		N[2] = 0.5 * xi * (xi + 1.0);
	}
	else if constexpr (Order == 3){
		N[0] = (-9.0/16.0) * (xi - (1.0/3.0)) * (xi - 1.0) * (xi + (1.0/3.0));
		N[1] = (27.0/16.0) * (xi - 1.0) * (xi - (1.0/3.0)) * (xi + 1.0);
		N[2] = (-27.0/16.0) * (xi + (1.0/3.0)) * (xi - 1.0) * (xi + 1.0);
		N[3] = (9.0/16.0) * (xi - (1.0/3.0)) * (xi + 1.0) * (xi + (1.0/3.0));
	}
	else {
		static_assert(Order <= 3, "Lagrange1D: Orders 1,2,3 implemented");
	}
}

// Implementation: evalFirstDerivative
template <int Order>
PDE_HOST PDE_DEVICE void Lagrange1D<Order>::evalFirstDerivative(Real xi, Real* N){
	if constexpr (Order == 1){
		N[0] = -0.5;
		N[1] = 0.5;
	}
	else if constexpr (Order == 2){
		N[0] = xi - 0.5;
		N[1] = -2.0 * xi;
		N[2] = xi + 0.5;
	}
	else if constexpr (Order == 3){
		N[0] = (-9.0/16.0) * (3.0*xi*xi - 2.0*xi - (1.0/9.0));
		N[1] = (27.0/16.0) * (3.0*xi*xi - (2.0/3.0)*xi - 1.0);
		N[2] = (-27.0/16.0) * (3.0*xi*xi + (2.0/3.0)*xi - 1.0);
		N[3] = (9.0/16.0) * (3.0*xi*xi + 2.0*xi - (1.0/9.0));
	}
	else {
		static_assert(Order <= 3, "Lagrange1D: Orders 1,2,3 implemented");
	}
}

// Implementation: evalSecondDerivative
template <int Order>
PDE_HOST PDE_DEVICE void Lagrange1D<Order>::evalSecondDerivative(Real xi, Real* N){
	if constexpr (Order == 1){
		N[0] = 0.0;
		N[1] = 0.0;
	}
	else if constexpr (Order == 2){
		N[0] = 1.0;
		N[1] = -2.0;
		N[2] = 1.0;
	}
	else if constexpr (Order == 3){
		N[0] = (-9.0/16.0) * (6.0*xi - 2.0);
		N[1] = (27.0/16.0) * (6.0*xi - (2.0/3.0));
		N[2] = (-27.0/16.0) * (6.0*xi + (2.0/3.0));
		N[3] = (9.0/16.0) * (6.0*xi + 2.0);
	}
	else {
		static_assert(Order <= 3, "Lagrange1D: Orders 1,2,3 implemented");
	}
}

} // namespace pdesolver::fem::basis
