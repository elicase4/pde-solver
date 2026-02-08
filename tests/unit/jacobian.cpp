#include <gtest/gtest.h>
#include <math.h>

#include "fem/basis/LagrangeQuad.hpp"
#include "fem/basis/LagrangeHex.hpp"
#include "fem/geometry/JacobianTransform.hpp"

using namespace pdesolver::fem::basis;
using namespace pdesolver::fem::geometry;

TEST(JacobianTransform, LagrangeQuad2Dx2DMainTransform){

	const Index pD = 2;
	const Index sD = 2;
	const Index npe = 4;

	Real nodeCoords[sD*npe] = {
		0, 0,
		2, 0,
		2, 3,
		0, 3
	};

	Real xi[pD] = {-0.2, -0.4};
	
	Real dNdxi[pD*npe];
	BilinearQuad::evalGradient(xi, dNdxi);

	Real J[sD*pD];
	JacobianTransform<sD,pD,npe>::computeJacobian(nodeCoords, dNdxi, J);
	EXPECT_NEAR(J[0], -xi[1], 1e-13);
	EXPECT_NEAR(J[1], -xi[0], 1e-13);
	EXPECT_NEAR(J[2], 0.0, 1e-13);
	EXPECT_NEAR(J[3], 1.5, 1e-13);

	Real g[pD*pD];
	JacobianTransform<sD,pD,npe>::computeMetric(J, g);
	EXPECT_NEAR(g[0], xi[1] * xi[1], 1e-13);
	EXPECT_NEAR(g[1], xi[1] * xi[0], 1e-13);
	EXPECT_NEAR(g[2], xi[1] * xi[0], 1e-13);
	EXPECT_NEAR(g[3], xi[0] * xi[0] + 2.25, 1e-13);

	Real measure = JacobianTransform<sD,pD,npe>::computeMeasure(g);
	EXPECT_NEAR(measure, sqrt((xi[1] * xi[1] * (xi[0] * xi[0] + 2.25) ) - (xi[1] * xi[1] * xi[0]* xi[0])), 1e-13);

}

TEST(JacobianTransform, LagrangeQuad2Dx2DGradientTransform){

	const Index pD = 2;
	const Index sD = 2;
	const Index npe = 4;

	Real nodeCoords[sD*npe] = {
		0, 0,
		2, 0,
		2, 3,
		0, 3
	};

	Real xi[pD] = {-0.2, -0.4};
	
	Real dNdxi[pD*npe];
	BilinearQuad::evalGradient(xi, dNdxi);

	Real J[sD*pD];
	JacobianTransform<sD,pD,npe>::computeJacobian(nodeCoords, dNdxi, J);
	
	Real g[pD*pD];
	JacobianTransform<sD,pD,npe>::computeMetric(J, g);

	Real dNdx[sD*npe];
	JacobianTransform<sD,pD,npe>::transformGradient(J, g, dNdxi, dNdx);

	Real sx=0, sy=0;
	for (Index a = 0; a < npe; ++a){
		sx += dNdx[a*sD   ];
		sy += dNdx[a*sD + 1];
	}

	EXPECT_NEAR(sx, 0.0, 1e-13);
	EXPECT_NEAR(sy, 0.0, 1e-13);

}

TEST(JacobianTransform, LagrangeQuad3Dx2DMainTransform){
	
	const Index pD = 2;
	const Index sD = 3;
	const Index npe = 4;

	Real nodeCoords[sD*npe] = {
		-2, -3,  4,
		 2, -3, -2,
		-2,  3,  3,
		 2,  3, -1
	};
	
	Real xi[pD] = {0.2, -0.1};
	
	Real dNdxi[pD*npe];
	BilinearQuad::evalGradient(xi, dNdxi);
	
	Real J[sD*pD];
	JacobianTransform<sD,pD,npe>::computeJacobian(nodeCoords, dNdxi, J);
	EXPECT_NEAR(J[0], 2.0, 1e-13);
	EXPECT_NEAR(J[1], 0.0, 1e-13);
	EXPECT_NEAR(J[2], 0.0, 1e-13);
	EXPECT_NEAR(J[3], 3.0, 1e-13);
	EXPECT_NEAR(J[4], -2.5+0.5*xi[1], 1e-13);
	EXPECT_NEAR(J[5], 0.5*xi[0], 1e-13);
	
	Real g[pD*pD];
	JacobianTransform<sD,pD,npe>::computeMetric(J, g);
	EXPECT_NEAR(g[0], 10.25 - 2.5*xi[1] + 0.25*xi[1]*xi[1], 1e-13);
	EXPECT_NEAR(g[1], -1.25*xi[0] + 0.25*xi[0]*xi[1], 1e-13);
	EXPECT_NEAR(g[2], -1.25*xi[0] + 0.25*xi[0]*xi[1], 1e-13);
	EXPECT_NEAR(g[3], 9.0 + 0.25*xi[0]*xi[0], 1e-13);

	Real measure = JacobianTransform<sD,pD,npe>::computeMeasure(g);
	EXPECT_NEAR(measure, sqrt(92.25 + xi[0]*xi[0] - 22.5*xi[1] + 2.25*xi[1]*xi[1]), 1e-13);

}

TEST(JacobianTransform, LagrangeQuad3Dx2DBoundaryNormal){
	
	const Index pD = 2;
	const Index sD = 3;
	const Index npe = 4;

	Real nodeCoords[sD*npe] = {
		-2, -3,  4,
		 2, -3, -2,
		-2,  3,  3,
		 2,  3, -1
	};
	
	Real xi[pD] = {-0.5, -1.0};

	Real dNdxi[pD*npe];
	BilinearQuad::evalGradient(xi, dNdxi);
	
	Real J[sD*pD];
	JacobianTransform<sD,pD,npe>::computeJacobian(nodeCoords, dNdxi, J);
	
	Index tangentID[pD];
	Real nCoeff = BilinearQuad::getFaceTopology(0, tangentID);
	EXPECT_NEAR(nCoeff, -1.0, 1e-13);
	EXPECT_NEAR(tangentID[0], 0, 1e-13);
	EXPECT_NEAR(tangentID[1], 1, 1e-13);
	
	Real normal[sD];
	JacobianTransform<sD,pD,npe>::computeNormal(J, tangentID, nCoeff, normal);

	EXPECT_NEAR(normal[0], -1.0 * (-1.25 * xi[0] + 0.25*xi[0]*xi[1])*2.0, 1e-13);
	EXPECT_NEAR(normal[1], (10.25 - 2.5*xi[1] + 0.25*xi[1]*xi[1])*3.0, 1e-13);
	EXPECT_NEAR(normal[2], -1.0 * (-1.25 * xi[0] + 0.25*xi[0]*xi[1])*(-2.5 + 0.5*xi[1]) + (10.25 - 2.5*xi[1] + 0.25*xi[1]*xi[1])*(0.5*xi[0]), 1e-13);
}
