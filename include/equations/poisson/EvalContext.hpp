#ifndef POISSON_EVALCONTEXT_HPP
#define POISSON_EVALCONTEXT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/EvalContext.hpp"

namespace pdesolver::fem::eval {

	template<typename Geometry, typename Basis, Int Dim, Int NodesPerElement>
	class PoissonEvalContext;

	template<typename Geometry, typename Basis, Int NodesPerElement>
	class PoissonEvalContext<Geometry, Basis, 2, NodesPerElement> {
	public:
		static constexpr Int Dimension = 2;
		static constexpr Int NumNodes = NodesPerElement;
		
		// node coordinates
		const Real* nodeCoords;
		
		// physical coordinate
		Real x[2];

		// time coordinate
		Real t;

		// quadrature
		Real xi[2];
		Real w;

		// geometry
		Real J[4];
		Real invJ[4];
		Real detJ;

		// basis values
		Real N[NodesPerElement];

		// ref gradients
		Real dNdxi[NodesPerElement];
		Real dNdeta[NodesPerElement];

		// physical gradients
		Real dNdx[NodesPerElement];
		Real dNdy[NodesPerElement];

		PDE_HOST PDE_DEVICE bindElement(const Real* coords, const Real time){
			nodeCoords = coords;
			t = time;
		}

		PDE_HOST PDE_DEVICE evaluate(const Real* xi_q, const Real weight){
			// set quad info
			xi[0] = xi_q[0];
			xi[1] = xi_q[1];
			w = weight;
			
			// evaluate basis
			Basis::evaluate(xi, N);
			Basis::evaluateGradient(xi, dNdxi, dNdeta);

			// geometry
			detJ = Geoemtry::computeMetric(nodeCoords, dNdxi, dNdeta, J);
			Geometry::invertMetric(J, detJ, invJ);
			Geometry::mapToPhysical(nodeCoords, N, x);

			// transforms
			Geometry::transformGradient(invJ, dNdxi, dNdeta, dNdx, dNdy);
		}

	}; // class PoissonEvalContext<Geometry, Basis, 2, NodesPerElement>
	
	template<typename Geometry, typename Basis, Int NodesPerElement>
	class PoissonEvalContext<Geometry, Basis, 3, NodesPerElement> {
	public:
		static constexpr Int Dimension = 3;
		static constexpr Int NumNodes = NodesPerElement;
		
		// node coordinates
		const Real* nodeCoords;
		
		// physical coordinate
		Real x[3];

		// time coordinate
		Real t;

		// quadrature
		Real xi[3];
		Real w;

		// geometry
		Real J[9];
		Real invJ[9];
		Real detJ;

		// basis values
		Real N[NodesPerElement];

		// ref gradients
		Real dNdxi[NodesPerElement];
		Real dNdeta[NodesPerElement];

		// physical gradients
		Real dNdx[NodesPerElement];
		Real dNdy[NodesPerElement];

		PDE_HOST PDE_DEVICE bindElement(const Real* coords, const Real time){
			nodeCoords = coords;
			t = time;
		}

		PDE_HOST PDE_DEVICE evaluate(const Real* xi_q, const Real weight){
			
			// set quad info
			xi[0] = xi_q[0];
			xi[1] = xi_q[1];
			xi[2] = xi_q[2];
			w = weight;
			
			// evaluate basis
			Basis::evaluate(xi, N);
			Basis::evaluateGradient(xi, dNdxi, dNdeta, dNdzeta);

			// geometry
			detJ = Geoemtry::computeMetric(nodeCoords, dNdxi, dNdeta, dNdzeta, J);
			Geometry::invertMetric(J, detJ, invJ);
			Geometry::mapToPhysical(nodeCoords, N, x);

			// transforms
			Geometry::transformGradient(invJ, dNdxi, dNdeta, dNdzeta, dNdx, dNdy, dNdz);
		}

	}; // class PoissonEvalContext<Geometry, Basis, 3, NodesPerElement>


} // namespace pdesolver::fem::eval

#endif
