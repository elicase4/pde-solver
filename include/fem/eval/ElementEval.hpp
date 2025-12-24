#ifndef PDESOLVER_ELEMENTEVAL_HPP
#define PDESOLVER_ELEMENTEVAL_HPP

#include "core/Types.hpp"
#include "core/CudaMacros.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<Int Dim, Int NodesPerElement>
			class ElementEval; // class ElementEval
			
			// 2D specialization
			template<Int NodesPerElement>
			class ElementEval<2, NodesPerElement> {
			public:
				
				// Geometry
				Real* nodeCoords;
				
				// Jacobian Transform
				Real J[2*2];
				Real invJ[2*2];
				Real detJ;
				Real detInvJ;
				
				// Basis
				Real* N;
				
				// Basis Gradient
				Real* dNdxi;
				Real* dNdeta;
				
				// Basis Hessian
				Real* d2Nd2xi;
				Real* d2Nd2eta;
				Real* d2Ndetadxi;

				// Basis Laplacian
				Real* lapN;
				
				// Quadrature
				Real* xi;
				Real* w;
				
				// Constructor
				HOST_DEVICE ElementEval() : nodeCoords(nullptr),
											J(nullptr), invJ(nullptr), detJ(0.0), detInvJ(0.0),
											N(nullptr),
											dNdxi(nullptr), dNdeta(nullptr),
											d2Nd2xi(nullptr), d2Nd2eta(nullptr), d2Ndetadxi(nullptr),
											lapN(nullptr),
											xi(nullptr), w(nullptr) {}
			}
			
			// 3D specialization
			template<Int NodesPerElement>
			class ElementEval<3, NodesPerElement> {
			public:
				
				// Geometry
				Real* nodeCoords;
				
				// Jacobian Transform
				Real J[3*3];
				Real invJ[3*3];
				Real detJ;
				Real detInvJ;
				
				// Basis
				Real* N;
				
				// Basis Gradient
				Real* dNdxi;
				Real* dNdeta;
				Real* dNdzeta;
				
				// Basis Hessian
				Real* d2Nd2xi;
				Real* d2Nd2eta;
				Real* d2Nd2zeta;
				Real* d2Ndetadxi;
				Real* d2Ndetadzeta;
				Real* d2Ndxidzeta;

				// Basis Laplacian
				Real* lapN;
				
				// Quadrature
				Real* xi;
				Real* w;
				
				// Constructor
				HOST_DEVICE ElementEval() : nodeCoords(nullptr),
											J(nullptr), invJ(nullptr), detJ(0.0), detInvJ(0.0),
											N(nullptr),
											dNdxi(nullptr), dNdeta(nullptr), dNdzeta(nullptr),
											d2Nd2xi(nullptr), d2Nd2eta(nullptr), d2Nd2zeta(nullptr), d2Ndetadxi(nullptr), d2Ndetadzeta(nullptr), d2Ndxidzeta(nullptr),
											lapN(nullptr),
											xi(nullptr), w(nullptr) {}
				}

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
