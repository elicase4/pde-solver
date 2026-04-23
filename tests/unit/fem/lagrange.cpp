#include <gtest/gtest.h>

#include "fem/basis/Lagrange1D.hpp"
#include "fem/basis/LagrangeQuad.hpp"
#include "fem/basis/LagrangeHex.hpp"

using namespace pdesolver::fem::basis;

// ======================================================
// Lagrange1D Tests
// ======================================================

TEST(Lagrange1D, PartitionOfUnity) {
	constexpr int numTestPoints = 5;
	Real testPoints[numTestPoints] = {-1.0, -0.5, 0.0, 0.5, 1.0};
	Real sum;

	for (int q = 0; q < numTestPoints; ++q){
		Real N[2];
		Lagrange1D<1>::eval(testPoints[q], N);
		sum = N[0] + N[1];
		EXPECT_NEAR(sum, 1.0, 1e-14) << "eval() p = " << 1 << " at xi = " << testPoints[q];
	}
	
	for (int q = 0; q < numTestPoints; ++q){
		Real N[3];
		Lagrange1D<2>::eval(testPoints[q], N);
		sum = N[0] + N[1] + N[2];
		EXPECT_NEAR(sum, 1.0, 1e-14) << "eval() p = " << 2 << " at xi = " << testPoints[q];
	}
	
	for (int q = 0; q < numTestPoints; ++q){
		Real N[4];
		Lagrange1D<3>::eval(testPoints[q], N);
		sum = N[0] + N[1] + N[2] + N[3];
		EXPECT_NEAR(sum, 1.0, 1e-14) << "eval() p = " << 3 << " at xi = " << testPoints[q];
	}
	
}

TEST(Lagrange1D, PartitionOfUnityDeriv) {
	constexpr int numTestPoints = 5;
	Real testPoints[numTestPoints] = {-1.0, -0.5, 0.0, 0.5, 1.0};
	Real sum;

	for (int q = 0; q < numTestPoints; ++q){
		Real N[2];
		Lagrange1D<1>::evalFirstDerivative(testPoints[q], N);
		sum = N[0] + N[1];
		EXPECT_NEAR(sum, 0.0, 1e-14) << "evalFirstDerivative() p = " << 1 << " at xi = " << testPoints[q];
		Lagrange1D<1>::evalSecondDerivative(testPoints[q], N);
		sum = N[0] + N[1];
		EXPECT_NEAR(sum, 0.0, 1e-14) << "evalSecondDerivative() p = " << 1 << " at xi = " << testPoints[q];
	}
	
	for (int q = 0; q < numTestPoints; ++q){
		Real N[3];
		Lagrange1D<2>::evalFirstDerivative(testPoints[q], N);
		sum = N[0] + N[1] + N[2];
		EXPECT_NEAR(sum, 0.0, 1e-14) << "evalFirstDerivative() p = " << 2 << " at xi = " << testPoints[q];
		Lagrange1D<2>::evalSecondDerivative(testPoints[q], N);
		sum = N[0] + N[1] + N[2];
		EXPECT_NEAR(sum, 0.0, 1e-14) << "evalSecondDerivative() p = " << 2 << " at xi = " << testPoints[q];
	}
	
	for (int q = 0; q < numTestPoints; ++q){
		Real N[4];
		Lagrange1D<3>::evalFirstDerivative(testPoints[q], N);
		sum = N[0] + N[1] + N[2] + N[3];
		EXPECT_NEAR(sum, 0.0, 1e-14) << "evalFirstDerivative() p = " << 3 << " at xi = " << testPoints[q];
		Lagrange1D<3>::evalSecondDerivative(testPoints[q], N);
		sum = N[0] + N[1] + N[2] + N[3];
		EXPECT_NEAR(sum, 0.0, 1e-14) << "evalSecondDerivative() p = " << 3 << " at xi = " << testPoints[q];
	}
	
}

TEST(Lagrange1D, KroneckerDeltaProperty) {
	
	Real nodes_linear[2] = {-1.0, 1.0};
	for (int i = 0; i < 2; ++i){
		Real N[2];
		Lagrange1D<1>::eval(nodes_linear[i], N);
		for (int j = 0; j < 2; ++j){
			Real expected = (i == j) ? 1.0 : 0.0;
			EXPECT_NEAR(N[j], expected, 1e-14) << "p = 1: N [" << j << "] at node " << i;
		}
	}

	Real nodes_quadratic[3] = {-1.0, 0.0, 1.0};
	for (int i = 0; i < 3; ++i){
		Real N[3];
		Lagrange1D<2>::eval(nodes_quadratic[i], N);
		for (int j = 0; j < 3; ++j){
			Real expected = (i == j) ? 1.0 : 0.0;
			EXPECT_NEAR(N[j], expected, 1e-14) << "p = 2: N [" << j << "] at node " << i;
		}
	}

	Real nodes_cubic[4] = {-1.0, (-1.0/3.0), (1.0/3.0), 1.0};
	for (int i = 0; i < 4; ++i){
		Real N[4];
		Lagrange1D<3>::eval(nodes_cubic[i], N);
		for (int j = 0; j < 4; ++j){
			Real expected = (i == j) ? 1.0 : 0.0;
			EXPECT_NEAR(N[j], expected, 1e-14) << "p = 3: N [" << j << "] at node " << i;
		}
	}

}

// ======================================================
// LagrangeQuad Tests
// ======================================================

TEST(LagrangeQuad, PartitionOfUnity) {
	constexpr int numTestPoints = 5;
	Real testPoints[numTestPoints][2] = {{-0.9,-0.8}, {0.0, -1.0}, {0.75, 0.6}, {0.5, -0.3}, {-0.7, 0.8}};
	Real sum;

	for (int q = 0; q < numTestPoints; ++q){
		Real N[4];
		BilinearQuad::eval(testPoints[q], N);
		sum = 0.0;
		for (int i = 0; i < 4; ++i) sum += N[i];
		EXPECT_NEAR(sum, 1.0, 1e-14) << "eval() p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
	}
	
	for (int q = 0; q < numTestPoints; ++q){
		Real N[9];
		BiquadraticQuad::eval(testPoints[q], N);
		sum = 0.0;
		for (int i = 0; i < 9; ++i) sum += N[i];
		EXPECT_NEAR(sum, 1.0, 1e-14) << "eval() p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
	}
	
	for (int q = 0; q < numTestPoints; ++q){
		Real N[16];
		BicubicQuad::eval(testPoints[q], N);
		sum = 0.0;
		for (int i = 0; i < 16; ++i) sum += N[i];
		EXPECT_NEAR(sum, 1.0, 1e-14) << "eval() p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
	}
	
}

TEST(LagrangeQuad, PartitionOfUnityDeriv) {
	constexpr int numTestPoints = 5;
	Real testPoints[numTestPoints][2] = {{-0.9,-0.8}, {0.0, -1.0}, {0.75, 0.6}, {0.5, -0.3}, {-0.7, 0.8}};

	for (int q = 0; q < numTestPoints; ++q){
		
		// Gradient
        Real dNdxi[8];
        BilinearQuad::evalGradient(testPoints[q], dNdxi);
        Real sum_dxi = 0.0;
        Real sum_deta = 0.0;
		for (int i = 0; i < 4; ++i){
			sum_dxi += dNdxi[2*i];
			sum_deta += dNdxi[2*i+1];
		}
        EXPECT_NEAR(sum_dxi, 0.0, 1e-14) << "evalGradient().xi p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_deta, 0.0, 1e-14) << "evalGradient().eta p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        
        // Hessian
        Real d2Nd2xi[12];
        BilinearQuad::evalHessian(testPoints[q], d2Nd2xi);
        Real sum_d2xi = 0.0;
        Real sum_d2eta = 0.0;
        Real sum_dxideta = 0.0;
		for (int i = 0; i < 4; ++i){ 
			sum_d2xi += d2Nd2xi[3*i];
			sum_d2eta += d2Nd2xi[3*i+1];
			sum_dxideta += d2Nd2xi[3*i+2];
		}
        EXPECT_NEAR(sum_d2xi, 0.0, 1e-14) << "evalHessian().xi2 p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_d2eta, 0.0, 1e-14) << "evalHessian().eta2 p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_dxideta, 0.0, 1e-14) << "evalHessian().xieta p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        
        // Laplacian
        Real lapN[4];
        BilinearQuad::evalLaplacian(testPoints[q], lapN);
        Real sum_lapN = 0.0;
		for (int i = 0; i < 4; ++i) sum_lapN += lapN[i];
        EXPECT_NEAR(sum_lapN, 0.0, 1e-14) << "evalLaplacian() p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
	}
	
	for (int q = 0; q < numTestPoints; ++q){
		
		// Gradient
        Real dNdxi[18];
        BiquadraticQuad::evalGradient(testPoints[q], dNdxi);
        Real sum_dxi = 0.0;
        Real sum_deta = 0.0;
		for (int i = 0; i < 9; ++i){
			sum_dxi += dNdxi[2*i];
			sum_deta += dNdxi[2*i+1];
		}
        EXPECT_NEAR(sum_dxi, 0.0, 1e-14) << "evalGradient().xi p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_deta, 0.0, 1e-14) << "evalGradient().eta p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        
        // Hessian
        Real d2Nd2xi[27];
        BiquadraticQuad::evalHessian(testPoints[q], d2Nd2xi);
        Real sum_d2xi = 0.0;
        Real sum_d2eta = 0.0;
        Real sum_dxideta = 0.0;
		for (int i = 0; i < 9; ++i){ 
			sum_d2xi += d2Nd2xi[3*i];
			sum_d2eta += d2Nd2xi[3*i+1];
			sum_dxideta += d2Nd2xi[3*i+2];
		}
        EXPECT_NEAR(sum_d2xi, 0.0, 1e-14) << "evalHessian().xi2 p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_d2eta, 0.0, 1e-14) << "evalHessian().eta2 p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_dxideta, 0.0, 1e-14) << "evalHessian().xieta p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        
        // Laplacian
        Real lapN[9];
        BiquadraticQuad::evalLaplacian(testPoints[q], lapN);
        Real sum_lapN = 0.0;
		for (int i = 0; i < 9; ++i) sum_lapN += lapN[i];
        EXPECT_NEAR(sum_lapN, 0.0, 1e-14) << "evalLaplacian() p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
	}

	for (int q = 0; q < numTestPoints; ++q){
		
		// Gradient
        Real dNdxi[32];
        BicubicQuad::evalGradient(testPoints[q], dNdxi);
        Real sum_dxi = 0.0;
        Real sum_deta = 0.0;
		for (int i = 0; i < 16; ++i){
			sum_dxi += dNdxi[2*i];
			sum_deta += dNdxi[2*i+1];
		}
        EXPECT_NEAR(sum_dxi, 0.0, 1e-14) << "evalGradient().xi p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_deta, 0.0, 1e-14) << "evalGradient().eta p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        
        // Hessian
        Real d2Nd2xi[48];
        BicubicQuad::evalHessian(testPoints[q], d2Nd2xi);
        Real sum_d2xi = 0.0;
        Real sum_d2eta = 0.0;
        Real sum_dxideta = 0.0;
		for (int i = 0; i < 16; ++i){ 
			sum_d2xi += d2Nd2xi[3*i];
			sum_d2eta += d2Nd2xi[3*i+1];
			sum_dxideta += d2Nd2xi[3*i+2];
		}
        EXPECT_NEAR(sum_d2xi, 0.0, 1e-14) << "evalHessian().xi2 p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_d2eta, 0.0, 1e-14) << "evalHessian().eta2 p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_dxideta, 0.0, 1e-14) << "evalHessian().xieta p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        
        // Laplacian
        Real lapN[16];
        BicubicQuad::evalLaplacian(testPoints[q], lapN);
        Real sum_lapN = 0.0;
		for (int i = 0; i < 16; ++i) sum_lapN += lapN[i];
        EXPECT_NEAR(sum_lapN, 0.0, 1e-14) << "evalLaplacian() p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
	}

}

// ======================================================
// LagrangeHex Tests
// ======================================================

TEST(LagrangeHex, PartitionOfUnity) {
	constexpr int numTestPoints = 5;
	Real testPoints[numTestPoints][3] = {{-0.9,-0.8,-0.6}, {0.0, -1.0, 0.25}, {0.75, 0.6, 0.9}, {0.5, -0.3, -0.6}, {-0.7, 0.8, 0.0}};
	Real sum;

	for (int q = 0; q < numTestPoints; ++q){
		Real N[8];
		TrilinearHex::eval(testPoints[q], N);
		sum = 0.0;
		for (int i = 0; i < 8; ++i) sum += N[i];
		EXPECT_NEAR(sum, 1.0, 1e-14) << "eval() p = (" << 1 << "," << 1 << ","<< 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][1] << ")";
	}
	
	for (int q = 0; q < numTestPoints; ++q){
		Real N[27];
		TriquadraticHex::eval(testPoints[q], N);
		sum = 0.0;
		for (int i = 0; i < 27; ++i) sum += N[i];
		EXPECT_NEAR(sum, 1.0, 1e-14) << "eval() p = (" << 2 << "," << 2 << ","<< 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][1] << ")";
	}
	
	for (int q = 0; q < numTestPoints; ++q){
		Real N[64];
		TricubicHex::eval(testPoints[q], N);
		sum = 0.0;
		for (int i = 0; i < 64; ++i) sum += N[i];
		EXPECT_NEAR(sum, 1.0, 1e-14) << "eval() p = (" << 3 << "," << 3 << ","<< 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][1] << ")";
	}
	
}

TEST(LagrangeHex, PartitionOfUnityDeriv) {
	constexpr int numTestPoints = 5;
	Real testPoints[numTestPoints][3] = {{-0.9,-0.8,-0.6}, {0.0, -1.0, 0.25}, {0.75, 0.6, 0.9}, {0.5, -0.3, -0.6}, {-0.7, 0.8, 0.0}};

	for (int q = 0; q < numTestPoints; ++q){
		
		// Gradient
        Real dNdxi[24];
        TrilinearHex::evalGradient(testPoints[q], dNdxi);
        Real sum_dxi = 0.0;
        Real sum_deta = 0.0;
        Real sum_dzeta = 0.0;
		for (int i = 0; i < 8; ++i){
			sum_dxi += dNdxi[3*i];
			sum_deta += dNdxi[3*i+1];
			sum_dzeta += dNdxi[3*i+2];
		}
        EXPECT_NEAR(sum_dxi, 0.0, 1e-14) << "evalGradient().xi p = (" << 1 << "," << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_deta, 0.0, 1e-14) << "evalGradient().eta p = (" << 1 << "," << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_dzeta, 0.0, 1e-14) << "evalGradient().zeta p = (" << 1 << "," << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        
        // Hessian
        Real d2Nd2xi[48];
        TrilinearHex::evalHessian(testPoints[q], d2Nd2xi);
        Real sum_d2xi = 0.0;
        Real sum_d2eta = 0.0;
        Real sum_d2zeta = 0.0;
        Real sum_detadxi = 0.0;
        Real sum_detadzeta = 0.0;
        Real sum_dxidzeta = 0.0;
		for (int i = 0; i < 8; ++i){ 
			sum_d2xi += d2Nd2xi[6*i];
			sum_d2eta += d2Nd2xi[6*i+1];
			sum_d2zeta += d2Nd2xi[6*i+2];
			sum_detadxi += d2Nd2xi[6*i+3];
			sum_detadzeta += d2Nd2xi[6*i+4];
			sum_dxidzeta += d2Nd2xi[6*i+5];
		}
        EXPECT_NEAR(sum_d2xi, 0.0, 1e-14) << "evalHessian().xi2 p = (" << 1 << "," << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_d2eta, 0.0, 1e-14) << "evalHessian().eta2 p = (" << 1 << "," << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_d2zeta, 0.0, 1e-14) << "evalHessian().zeta2 p = (" << 1 << "," << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_detadxi, 0.0, 1e-14) << "evalHessian().etaxi p = (" << 1 << "," << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_detadzeta, 0.0, 1e-14) << "evalHessian().etazeta p = (" << 1 << "," << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_dxidzeta, 0.0, 1e-14) << "evalHessian().xizeta p = (" << 1 << "," << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        
        // Laplacian
        Real lapN[8];
        TrilinearHex::evalLaplacian(testPoints[q], lapN);
        Real sum_lapN = 0.0;
		for (int i = 0; i < 8; ++i) sum_lapN += lapN[i];
        EXPECT_NEAR(sum_d2xi, 0.0, 1e-14) << "evalLaplacian() p = (" << 1 << "," << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
	}

	for (int q = 0; q < numTestPoints; ++q){
		
		// Gradient
        Real dNdxi[81];
        TriquadraticHex::evalGradient(testPoints[q], dNdxi);
        Real sum_dxi = 0.0;
        Real sum_deta = 0.0;
        Real sum_dzeta = 0.0;
		for (int i = 0; i < 27; ++i){
			sum_dxi += dNdxi[3*i];
			sum_deta += dNdxi[3*i+1];
			sum_deta += dNdxi[3*i+2];
		}
        EXPECT_NEAR(sum_dxi, 0.0, 1e-14) << "evalGradient().xi p = (" << 2 << "," << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_deta, 0.0, 1e-14) << "evalGradient().eta p = (" << 2 << "," << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_dzeta, 0.0, 1e-14) << "evalGradient().zeta p = (" << 2 << "," << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        
        // Hessian
        Real d2Nd2xi[162];
        TriquadraticHex::evalHessian(testPoints[q], d2Nd2xi);
        Real sum_d2xi = 0.0;
        Real sum_d2eta = 0.0;
        Real sum_d2zeta = 0.0;
        Real sum_detadxi = 0.0;
        Real sum_detadzeta = 0.0;
        Real sum_dxidzeta = 0.0;
		for (int i = 0; i < 27; ++i){ 
			sum_d2xi += d2Nd2xi[6*i];
			sum_d2eta += d2Nd2xi[6*i+1];
			sum_d2zeta += d2Nd2xi[6*i+2];
			sum_detadxi += d2Nd2xi[6*i+3];
			sum_detadzeta += d2Nd2xi[6*i+4];
			sum_dxidzeta += d2Nd2xi[6*i+5];
		}
        EXPECT_NEAR(sum_d2xi, 0.0, 1e-14) << "evalHessian().xi2 p = (" << 2 << "," << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_d2eta, 0.0, 1e-14) << "evalHessian().eta2 p = (" << 2 << "," << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_d2zeta, 0.0, 1e-14) << "evalHessian().zeta2 p = (" << 2 << "," << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_detadxi, 0.0, 1e-14) << "evalHessian().etaxi p = (" << 2 << "," << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_detadzeta, 0.0, 1e-14) << "evalHessian().etazeta p = (" << 2 << "," << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_dxidzeta, 0.0, 1e-14) << "evalHessian().xizeta p = (" << 2 << "," << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        
        // Laplacian
        Real lapN[27];
        TriquadraticHex::evalLaplacian(testPoints[q], lapN);
        Real sum_lapN = 0.0;
		for (int i = 0; i < 27; ++i) sum_lapN += lapN[i];
        EXPECT_NEAR(sum_d2xi, 0.0, 1e-14) << "evalLaplacian() p = (" << 2 << "," << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
	}

	for (int q = 0; q < numTestPoints; ++q){
		
		// Gradient
        Real dNdxi[192];
        TricubicHex::evalGradient(testPoints[q], dNdxi);
        Real sum_dxi = 0.0;
        Real sum_deta = 0.0;
        Real sum_dzeta = 0.0;
		for (int i = 0; i < 64; ++i){
			sum_dxi += dNdxi[3*i];
			sum_deta += dNdxi[3*i+1];
			sum_deta += dNdxi[3*i+2];
		}
        EXPECT_NEAR(sum_dxi, 0.0, 1e-14) << "evalGradient().xi p = (" << 3 << "," << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_deta, 0.0, 1e-14) << "evalGradient().eta p = (" << 3 << "," << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_dzeta, 0.0, 1e-14) << "evalGradient().zeta p = (" << 3 << "," << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        
        // Hessian
        Real d2Nd2xi[384];
        TricubicHex::evalHessian(testPoints[q], d2Nd2xi);
        Real sum_d2xi = 0.0;
        Real sum_d2eta = 0.0;
        Real sum_d2zeta = 0.0;
        Real sum_detadxi = 0.0;
        Real sum_detadzeta = 0.0;
        Real sum_dxidzeta = 0.0;
		for (int i = 0; i < 64; ++i){ 
			sum_d2xi += d2Nd2xi[6*i];
			sum_d2eta += d2Nd2xi[6*i+1];
			sum_d2zeta += d2Nd2xi[6*i+2];
			sum_detadxi += d2Nd2xi[6*i+3];
			sum_detadzeta += d2Nd2xi[6*i+4];
			sum_dxidzeta += d2Nd2xi[6*i+5];
		}
        EXPECT_NEAR(sum_d2xi, 0.0, 1e-14) << "evalHessian().xi2 p = (" << 3 << "," << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_d2eta, 0.0, 1e-14) << "evalHessian().eta2 p = (" << 3 << "," << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_d2zeta, 0.0, 1e-14) << "evalHessian().zeta2 p = (" << 3 << "," << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_detadxi, 0.0, 1e-14) << "evalHessian().etaxi p = (" << 3 << "," << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_detadzeta, 0.0, 1e-14) << "evalHessian().etazeta p = (" << 3 << "," << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        EXPECT_NEAR(sum_dxidzeta, 0.0, 1e-14) << "evalHessian().xizeta p = (" << 3 << "," << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
        
        // Laplacian
        Real lapN[64];
        TricubicHex::evalLaplacian(testPoints[q], lapN);
        Real sum_lapN = 0.0;
		for (int i = 0; i < 64; ++i) sum_lapN += lapN[i];
        EXPECT_NEAR(sum_d2xi, 0.0, 1e-14) << "evalLaplacian() p = (" << 3 << "," << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << "," << testPoints[q][2] << ")";
	}

}
