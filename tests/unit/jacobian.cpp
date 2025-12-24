#include <gtest/gtest.h>

#include "fem/basis/LagrangeQuad.hpp"
#include "fem/basis/LagrangeHex.hpp"
#include "fem/geometry/JacobianTransform.hpp"

using namespace pdesolver::fem::basis;
using namespace pdesolver::fem::geometry;

// ============================================
// 2D Jacobian Tests - Reference Domain [-1,1]^2
// ============================================

TEST(JacobianTransform2D, BilinearQuad) {
    // Bilinear quad with nodes at corners of reference domain [-1,1]^2
    // Mapped to physical: [0,2] x [0,3]
    Real nodeCoords[8] = {
        0.0, 0.0,  // Node 0: (-1,-1)
        2.0, 0.0,  // Node 1: ( 1,-1)
        0.0, 3.0,  // Node 2: (-1, 1)
        2.0, 3.0   // Node 3: ( 1, 1)
    };
    
    // Test at three points in reference domain
    Real testPoints[3][2] = {
        {-0.5, -0.5},
        { 0.0,  0.3},
        { 0.7,  0.2}
    };
    
    for (int q = 0; q < 3; ++q) {
        Real xi[2] = {testPoints[q][0], testPoints[q][1]};
        
        // Evaluate basis gradients
        Real dNdxi[4], dNdeta[4];
        BilinearQuad::evalGradient(xi, dNdxi, dNdeta);
        
        // Compute forward Jacobian
        Real J[4];
        Real detJ = JacobianTransform<2, 4>::computeForward(nodeCoords, dNdxi, dNdeta, J);
        
        // Check forward jacobian mapping
        EXPECT_NEAR(J[0], 1.0, 1e-13) << "Point " << q << ": dx/dxi";
        EXPECT_NEAR(J[1], 0.0, 1e-13) << "Point " << q << ": dx/deta";
        EXPECT_NEAR(J[2], 0.0, 1e-13) << "Point " << q << ": dy/dxi";
        EXPECT_NEAR(J[3], 1.5, 1e-13) << "Point " << q << ": dy/deta";
        
        // Determinant should be constant: 1.0 * 1.5 = 1.5
        EXPECT_NEAR(detJ, 1.5, 1e-13) << "Point " << q;
        
        // Invert Jacobian
        Real invJ[4];
        JacobianTransform<2, 4>::invertForward(J, detJ, invJ);
        
        // Check J * invJ = I
        Real I00 = J[0]*invJ[0] + J[1]*invJ[2];
        Real I01 = J[0]*invJ[1] + J[1]*invJ[3];
        Real I10 = J[2]*invJ[0] + J[3]*invJ[2];
        Real I11 = J[2]*invJ[1] + J[3]*invJ[3];
        
        EXPECT_NEAR(I00, 1.0, 1e-12) << "Point " << q;
        EXPECT_NEAR(I01, 0.0, 1e-12) << "Point " << q;
        EXPECT_NEAR(I10, 0.0, 1e-12) << "Point " << q;
        EXPECT_NEAR(I11, 1.0, 1e-12) << "Point " << q;
    }
}

TEST(JacobianTransform2D, PhysicalGradientTransform) {
    // Rectangle in physical space: [1,4] x [2,5]
    Real nodeCoords[8] = {
        1.0, 2.0,  // (-1,-1)
        4.0, 2.0,  // ( 1,-1)
        1.0, 5.0,  // (-1, 1)
        4.0, 5.0   // ( 1, 1)
    };
    
    Real testPoints[3][2] = {
        {-0.3, 0.5},
        { 0.0, 0.0},
        { 0.6, -0.4}
    };
    
    for (int q = 0; q < 3; ++q) {
        Real xi[2] = {testPoints[q][0], testPoints[q][1]};
        
        // Evaluate basis and gradients in reference space
        Real N[4];
        Real dNdxi[4], dNdeta[4];
        BilinearQuad::eval(xi, N);
        BilinearQuad::evalGradient(xi, dNdxi, dNdeta);
        
        // Compute Jacobian
        Real J[4];
        Real detJ = JacobianTransform<2, 4>::computeForward(nodeCoords, dNdxi, dNdeta, J);
        
        Real invJ[4];
        JacobianTransform<2, 4>::invertForward(J, detJ, invJ);
        
        // Transform gradients to physical space
        Real dNdx[4], dNdy[4];
        JacobianTransform<2, 4>::transformGradient(invJ, dNdxi, dNdeta, dNdx, dNdy);
        
        // Verify partition of unity: sum(dN/dx) = 0, sum(dN/dy) = 0
        Real sum_dNdx = 0.0, sum_dNdy = 0.0;
        for (int a = 0; a < 4; ++a) {
            sum_dNdx += dNdx[a];
            sum_dNdy += dNdy[a];
        }
        EXPECT_NEAR(sum_dNdx, 0.0, 1e-13) << "Point " << q;
        EXPECT_NEAR(sum_dNdy, 0.0, 1e-13) << "Point " << q;
        
        // Verify physical gradients are scaled correctly
        // J = [1.5  0 ], invJ = [2/3    0  ]
        //     [0   1.5]         [0    2/3]
        for (int a = 0; a < 4; ++a) {
            EXPECT_NEAR(dNdx[a], (2.0/3.0) * dNdxi[a], 1e-12) 
                << "Point " << q << ", node " << a;
            EXPECT_NEAR(dNdy[a], (2.0/3.0) * dNdeta[a], 1e-12) 
                << "Point " << q << ", node " << a;
        }
    }
}

// ============================================
// 3D Jacobian Tests - Reference Domain [-1,1]^3
// ============================================

TEST(JacobianTransform3D, TrilinearHex) {
    // Trilinear hex with nodes at corners of reference domain [-1,1]^3
    // Node ordering: standard hex ordering
    // Mapped to physical: [0,2] x [0,3] x [0,4]
    Real nodeCoords[24] = {
        0.0, 0.0, 0.0,  // Node 0: (-1,-1,-1)
        2.0, 0.0, 0.0,  // Node 1: ( 1,-1,-1)
        0.0, 3.0, 0.0,  // Node 2: (-1, 1,-1)
        2.0, 3.0, 0.0,  // Node 3: ( 1, 1,-1)
        0.0, 0.0, 4.0,  // Node 4: (-1,-1, 1)
        2.0, 0.0, 4.0,  // Node 5: ( 1,-1, 1)
        0.0, 3.0, 4.0,  // Node 6: (-1, 1, 1)
        2.0, 3.0, 4.0   // Node 7: ( 1, 1, 1)
    };
    
    Real testPoints[3][3] = {
        {-0.4, -0.3, 0.2},
        { 0.0,  0.0, 0.0},
        { 0.5,  0.7, -0.6}
    };
    
    for (int q = 0; q < 3; ++q) {
        Real xi[3] = {testPoints[q][0], testPoints[q][1], testPoints[q][2]};
        
        // Evaluate basis gradients
        Real dNdxi[8], dNdeta[8], dNdzeta[8];
        TrilinearHex::evalGradient(xi, dNdxi, dNdeta, dNdzeta);
        
        // Compute forward Jacobian
        Real J[9];
        Real detJ = JacobianTransform<3, 8>::computeForward(nodeCoords, dNdxi, dNdeta, dNdzeta, J);
        
        // For affine mapping, Jacobian should be constant
        // J = [1.0  0    0  ]
        //     [0    1.5  0  ]
        //     [0    0    2.0]
        EXPECT_NEAR(J[0], 1.0, 1e-13) << "Point " << q;
        EXPECT_NEAR(J[1], 0.0, 1e-13) << "Point " << q;
        EXPECT_NEAR(J[2], 0.0, 1e-13) << "Point " << q;
        EXPECT_NEAR(J[3], 0.0, 1e-13) << "Point " << q;
        EXPECT_NEAR(J[4], 1.5, 1e-13) << "Point " << q;
        EXPECT_NEAR(J[5], 0.0, 1e-13) << "Point " << q;
        EXPECT_NEAR(J[6], 0.0, 1e-13) << "Point " << q;
        EXPECT_NEAR(J[7], 0.0, 1e-13) << "Point " << q;
        EXPECT_NEAR(J[8], 2.0, 1e-13) << "Point " << q;
        
        // Determinant: 1.0 * 1.5 * 2.0 = 3.0
        EXPECT_NEAR(detJ, 3.0, 1e-13) << "Point " << q;
        
        // Invert Jacobian
        Real invJ[9];
        JacobianTransform<3, 8>::invertForward(J, detJ, invJ);
        
        // Check J * invJ = I (sample diagonal)
        Real I00 = J[0]*invJ[0] + J[1]*invJ[3] + J[2]*invJ[6];
        Real I11 = J[3]*invJ[1] + J[4]*invJ[4] + J[5]*invJ[7];
        Real I22 = J[6]*invJ[2] + J[7]*invJ[5] + J[8]*invJ[8];
        Real I01 = J[0]*invJ[1] + J[1]*invJ[4] + J[2]*invJ[7];
        
        EXPECT_NEAR(I00, 1.0, 1e-12) << "Point " << q;
        EXPECT_NEAR(I11, 1.0, 1e-12) << "Point " << q;
        EXPECT_NEAR(I22, 1.0, 1e-12) << "Point " << q;
        EXPECT_NEAR(I01, 0.0, 1e-12) << "Point " << q;
    }
}

TEST(JacobianTransform3D, PhysicalGradientTransform) {
    // Box: [1,3] x [2,5] x [0,6]
    Real nodeCoords[24] = {
        1.0, 2.0, 0.0,  // (-1,-1,-1)
        3.0, 2.0, 0.0,  // ( 1,-1,-1)
        1.0, 5.0, 0.0,  // (-1, 1,-1)
        3.0, 5.0, 0.0,  // ( 1, 1,-1)
        1.0, 2.0, 6.0,  // (-1,-1, 1)
        3.0, 2.0, 6.0,  // ( 1,-1, 1)
        1.0, 5.0, 6.0,  // (-1, 1, 1)
        3.0, 5.0, 6.0   // ( 1, 1, 1)
    };
    
    Real testPoints[3][3] = {
        {-0.2, 0.3, -0.5},
        { 0.0, 0.0,  0.0},
        { 0.4, -0.6, 0.8}
    };
    
    for (int q = 0; q < 3; ++q) {
        Real xi[3] = {testPoints[q][0], testPoints[q][1], testPoints[q][2]};
        
        // Evaluate basis and gradients
        Real dNdxi[8], dNdeta[8], dNdzeta[8];
        TrilinearHex::evalGradient(xi, dNdxi, dNdeta, dNdzeta);
        
        // Compute Jacobian
        Real J[9];
        Real detJ = JacobianTransform<3, 8>::computeForward(nodeCoords, dNdxi, dNdeta, dNdzeta, J);
        
        Real invJ[9];
        JacobianTransform<3, 8>::invertForward(J, detJ, invJ);
        
        // Transform gradients
        Real dNdx[8], dNdy[8], dNdz[8];
        JacobianTransform<3, 8>::transformGradient(invJ, dNdxi, dNdeta, dNdzeta, dNdx, dNdy, dNdz);
        
        // Verify partition of unity
        Real sum_dNdx = 0.0, sum_dNdy = 0.0, sum_dNdz = 0.0;
        for (int a = 0; a < 8; ++a) {
            sum_dNdx += dNdx[a];
            sum_dNdy += dNdy[a];
            sum_dNdz += dNdz[a];
        }
        EXPECT_NEAR(sum_dNdx, 0.0, 1e-12) << "Point " << q;
        EXPECT_NEAR(sum_dNdy, 0.0, 1e-12) << "Point " << q;
        EXPECT_NEAR(sum_dNdz, 0.0, 1e-12) << "Point " << q;
    }
}
