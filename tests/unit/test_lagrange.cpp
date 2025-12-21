#include <gtest/gtest.h>
#include <cmath>
#include "fem/Lagrange1D.hpp"
#include "fem/LagrangeQuad.hpp"
#include "fem/LagrangeHex.hpp"

using namespace pdesolver::fem;

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
        Real dNdxi[4], dNdtheta[4];
        BilinearQuad::evalGradient(testPoints[q], dNdxi, dNdtheta);
        Real sum_dxi = 0.0;
		for (int i = 0; i < 4; ++i) sum_dxi += dNdxi[i];
        Real sum_dtheta = 0.0;
		for (int i = 0; i < 4; ++i) sum_dtheta += dNdtheta[i];
        EXPECT_NEAR(sum_dxi, 0.0, 1e-14) << "evalGradient().xi p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_dtheta, 0.0, 1e-14) << "evalGradient().theta p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        
        // Hessian
        Real d2Nd2xi[4], d2Nd2theta[4], d2Ndxidtheta[4];
        BilinearQuad::evalHessian(testPoints[q], d2Nd2xi, d2Nd2theta, d2Ndxidtheta);
        Real sum_d2xi = 0.0;
		for (int i = 0; i < 4; ++i) sum_d2xi += d2Nd2xi[i];
        Real sum_d2theta = 0.0;
		for (int i = 0; i < 4; ++i) sum_d2theta += d2Nd2theta[i];
        Real sum_dxidtheta = 0.0;
		for (int i = 0; i < 4; ++i) sum_dxidtheta += d2Ndxidtheta[i];
        EXPECT_NEAR(sum_d2xi, 0.0, 1e-14) << "evalHessian().xi2 p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_d2theta, 0.0, 1e-14) << "evalHessian().theta2 p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_dxidtheta, 0.0, 1e-14) << "evalHessian().xitheta p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        
        // Laplacian
        Real lapN[4];
        BilinearQuad::evalLaplacian(testPoints[q], lapN);
        Real sum_lapN = 0.0;
		for (int i = 0; i < 4; ++i) sum_lapN += lapN[i];
        EXPECT_NEAR(sum_lapN, 0.0, 1e-14) << "evalLaplacian() p = (" << 1 << "," << 1 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
	}
	
	for (int q = 0; q < numTestPoints; ++q){
		
		// Gradient
        Real dNdxi[9], dNdtheta[9];
        BiquadraticQuad::evalGradient(testPoints[q], dNdxi, dNdtheta);
        Real sum_dxi = 0.0;
		for (int i = 0; i < 9; ++i) sum_dxi += dNdxi[i];
        Real sum_dtheta = 0.0;
		for (int i = 0; i < 9; ++i) sum_dtheta += dNdtheta[i];
        EXPECT_NEAR(sum_dxi, 0.0, 1e-14) << "evalGradient().xi p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_dtheta, 0.0, 1e-14) << "evalGradient().theta p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        
        // Hessian
        Real d2Nd2xi[9], d2Nd2theta[9], d2Ndxidtheta[9];
        BiquadraticQuad::evalHessian(testPoints[q], d2Nd2xi, d2Nd2theta, d2Ndxidtheta);
        Real sum_d2xi = 0.0;
		for (int i = 0; i < 9; ++i) sum_d2xi += d2Nd2xi[i];
        Real sum_d2theta = 0.0;
		for (int i = 0; i < 9; ++i) sum_d2theta += d2Nd2theta[i];
        Real sum_dxidtheta = 0.0;
		for (int i = 0; i < 9; ++i) sum_dxidtheta += d2Ndxidtheta[i];
        EXPECT_NEAR(sum_d2xi, 0.0, 1e-14) << "evalHessian().xi2 p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_d2theta, 0.0, 1e-14) << "evalHessian().theta2 p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_dxidtheta, 0.0, 1e-14) << "evalHessian().xitheta p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        
        // Laplacian
        Real lapN[9];
        BiquadraticQuad::evalLaplacian(testPoints[q], lapN);
        Real sum_lapN = 0.0;
		for (int i = 0; i < 9; ++i) sum_lapN += lapN[i];
        EXPECT_NEAR(sum_lapN, 0.0, 1e-14) << "evalLaplacian() p = (" << 2 << "," << 2 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
	}

	for (int q = 0; q < numTestPoints; ++q){
		
		// Gradient
        Real dNdxi[16], dNdtheta[16];
        BicubicQuad::evalGradient(testPoints[q], dNdxi, dNdtheta);
        Real sum_dxi = 0.0;
		for (int i = 0; i < 16; ++i) sum_dxi += dNdxi[i];
        Real sum_dtheta = 0.0;
		for (int i = 0; i < 16; ++i) sum_dtheta += dNdtheta[i];
        EXPECT_NEAR(sum_dxi, 0.0, 1e-14) << "evalGradient().xi p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_dtheta, 0.0, 1e-14) << "evalGradient().theta p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        
        // Hessian
        Real d2Nd2xi[16], d2Nd2theta[16], d2Ndxidtheta[16];
        BicubicQuad::evalHessian(testPoints[q], d2Nd2xi, d2Nd2theta, d2Ndxidtheta);
        Real sum_d2xi = 0.0;
		for (int i = 0; i < 16; ++i) sum_d2xi += d2Nd2xi[i];
        Real sum_d2theta = 0.0;
		for (int i = 0; i < 16; ++i) sum_d2theta += d2Nd2theta[i];
        Real sum_dxidtheta = 0.0;
		for (int i = 0; i < 16; ++i) sum_dxidtheta += d2Ndxidtheta[i];
        EXPECT_NEAR(sum_d2xi, 0.0, 1e-14) << "evalHessian().xi2 p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_d2theta, 0.0, 1e-14) << "evalHessian().theta2 p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        EXPECT_NEAR(sum_dxidtheta, 0.0, 1e-14) << "evalHessian().xitheta p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
        
        // Laplacian
        Real lapN[16];
        BicubicQuad::evalLaplacian(testPoints[q], lapN);
        Real sum_lapN = 0.0;
		for (int i = 0; i < 16; ++i) sum_lapN += lapN[i];
        EXPECT_NEAR(sum_lapN, 0.0, 1e-14) << "evalLaplacian() p = (" << 3 << "," << 3 << ") at xi = (" << testPoints[q][0] << "," << testPoints[q][1] << ")";
	}

}

