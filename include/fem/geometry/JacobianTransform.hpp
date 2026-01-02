#ifndef PDESOLVER_JACOBIANTRANSFORM_HPP
#define PDESOLVER_JACOBIANTRANSFORM_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace geometry {
			
			template<Int Dim, Int NodesPerElement>
			class JacobianTransform; // class JacobianTransform
			
			// 2D specialization
			template<Int NodesPerElement>
			class JacobianTransform<2, NodesPerElement> {
			public:
				
				// base transform
				PDE_HOST PDE_DEVICE static void computeForward(const Real* nodeCoords, const Real* N, Real* x);

				// jacobian matrix
				PDE_HOST PDE_DEVICE static Real computeJacobian(const Real* nodeCoords, const Real* dNdxi, const Real* dNdeta, Real* J);
				PDE_HOST PDE_DEVICE static Real invertJacobian(const Real* J, const Real detJ, Real* invJ);
				
				// operator transforms
				PDE_HOST PDE_DEVICE static void transformGradient(const Real* invJ, const Real* dNdxi, const Real* dNdeta, Real* dNdx, Real* dNdy);
				// TODO: add laplacian, hessian
				
				// normal computation
				PDE_HOST PDE_DEVICE static void computeNormal(const Real* J, const Index* tangentID, const Real nCoeff, Real* n);
			};
			
			// 3D specialization
			template<Int NodesPerElement>
			class JacobianTransform<3, NodesPerElement> {
			public:
				
				// base transform
				PDE_HOST PDE_DEVICE static void computeForward(const Real* nodeCoords, const Real* N, Real* x);

				// jacobian matrix
				PDE_HOST PDE_DEVICE static Real computeJacobian(const Real* nodeCoords, const Real* dNdxi, const Real* dNdeta, const Real* dNdzeta, Real* J);
				PDE_HOST PDE_DEVICE static Real invertJacobian(const Real* J, const Real detJ, Real* invJ);
				
				// operator transforms
				PDE_HOST PDE_DEVICE static void transformGradient(const Real* invJ, const Real* dNdxi, const Real* dNdeta, const Real* dNdzeta, Real* dNdx, Real* dNdy, Real* dNdz);
				// TODO: add laplacian, hessian
				
				// normal computation
				PDE_HOST PDE_DEVICE static void computeNormal(const Real* J, const Index* tangentID, const Real nCoeff, Real* n);
			};
			
			#include "JacobianTransform2D.tpp"
			#include "JacobianTransform3D.tpp"

		} // namespace geometry
	} // namespace fem
} // namespace pdesolver

#endif
