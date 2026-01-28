#ifndef PDESOLVER_ELEMENTEVAL_HPP
#define PDESOLVER_ELEMENTEVAL_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace eval {

			template<Int Dim, Int NodesPerElement>
			class EvalContext; // class EvalContext
			
			// 2D specialization
			template<Int NodesPerElement>
			class EvalContext<2, NodesPerElement> {
			public:
		
				// Parameters
				static constexpr Int Dimension = Dim;
				static constexpr Int NumNodes = NodesPerElement;
				
				// Geometry
				Real nodeCoords[NodesPerElement];
				Real* normal;

				// Time
				Real t;

				// Transform
				Real J[2*2];
				Real invJ[2*2];
				Real detJ;
				Real detInvJ;
				
				// Basis
				Real N[NodesPerElement];
				
				// Basis Gradient
				Real dNdxi[NodesPerElement];
				Real dNdeta[NodesPerElement];
				
				// Basis Divergence
				Real divN[NodesPerElement];
				
				// Basis Hessian
				Real d2Nd2xi[NodesPerElement];
				Real d2Nd2eta[NodesPerElement];
				Real d2Ndetadxi[NodesPerElement];

				// Basis Laplacian
				Real lapN[NodesPerElement];
				
				// Quadrature
				Real xi[2];
				Real w;
				
				// Constructor
				PDE_HOST PDE_DEVICE EvalContext() : nodeCoords(nullptr), normal(nullptr),
													J(nullptr), invJ(nullptr), detJ(0.0), detInvJ(0.0),
													N(nullptr),
													dNdxi(nullptr), dNdeta(nullptr),
													divN(nullptr),
													d2Nd2xi(nullptr), d2Nd2eta(nullptr), d2Ndetadxi(nullptr),
													lapN(nullptr),
													xi(nullptr), w(0.0) {}
			}
			
			// 3D specialization
			template<Int NodesPerElement>
			class EvalContext<3, NodesPerElement> {
			public:
		
				// Parameters
				static constexpr Int Dimension = Dim;
				static constexpr Int NumNodes = NodesPerElement;
				
				// Geometry
				Real nodeCoords[NodesPerElement];
				Real* normal;
				
				// Time
				Real t;
				
				// Transform
				Real J[3*3];
				Real invJ[3*3];
				Real detJ;
				Real detInvJ;
				
				// Basis
				Real N[NodesPerElement];
				
				// Basis Gradient
				Real dNdxi[NodesPerElement];
				Real dNdeta[NodesPerElement];
				Real dNdzeta[NodesPerElement];
				
				// Basis Divergence
				Real divN[NodesPerElement];
				
				// Basis Hessian
				Real d2Nd2xi[NodesPerElement];
				Real d2Nd2eta[NodesPerElement];
				Real d2Nd2zeta[NodesPerElement];
				Real d2Ndetadxi[NodesPerElement];
				Real d2Ndetadzeta[NodesPerElement];
				Real d2Ndxidzeta[NodesPerElement];

				// Basis Laplacian
				Real lapN[NodesPerElement];
				
				// Quadrature
				Real xi[3];
				Real w;
				
				// Constructor
				PDE_HOST PDE_DEVICE EvalContext() : nodeCoords(nullptr), normal(nullptr),
													J(nullptr), invJ(nullptr), detJ(0.0), detInvJ(0.0),
													N(nullptr),
													dNdxi(nullptr), dNdeta(nullptr), dNdzeta(nullptr),
													divN(nullptr),
													d2Nd2xi(nullptr), d2Nd2eta(nullptr), d2Nd2zeta(nullptr), d2Ndetadxi(nullptr), d2Ndetadzeta(nullptr), d2Ndxidzeta(nullptr),
													lapN(nullptr),
													xi(nullptr), w(0.0) {}
				}

		} // namespace eval
	} // namespace fem
} // namespace pdesolver

#endif
