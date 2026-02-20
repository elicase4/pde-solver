#ifndef PDESOLVER_JACOBIANTRANSFORM_HPP
#define PDESOLVER_JACOBIANTRANSFORM_HPP

#include <math.h>

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace geometry {
			
			template<Int SpatialDimension, Int ParametricDimension, Int NodesPerElement>
			class JacobianTransform {
			public:
				PDE_HOST PDE_DEVICE static void mapToPhysical(const Real* nodeCoords, const Real* N, Real* x);
				
				PDE_HOST PDE_DEVICE static void computeJacobian(const Real* nodeCoords, const Real* dNdxi, Real* J);

				PDE_HOST PDE_DEVICE static Real computeMeasure(const Real* g);

				PDE_HOST PDE_DEVICE static void computeMetric(const Real* J, Real* g);

				PDE_HOST PDE_DEVICE static void transformGradient(const Real* J, const Real* g, const Real* dNdxi, Real* dNdx);

				PDE_HOST PDE_DEVICE static void computeNormal(const Real* J, const Index* tangentID, const Real nCoeff, Real* n);

				Index SpatialDim = SpatialDimension;
				Index ParametricDim = ParametricDimension;
				Index NumNodes = NodesPerElement;
			
			private:
				PDE_HOST PDE_DEVICE static void invertMatrix(const Real* A, const Real detA, Real* invA);
				
				PDE_HOST PDE_DEVICE static Real computeMatrixDeterminant(const Real* A);

			}; // class JacobianTransform
			
		} // namespace geometry
	} // namespace fem
} // namespace pdesolver

#include "JacobianTransform.tpp"

#endif
