#ifndef POISSON_EVALCONTEXT_HPP
#define POISSON_EVALCONTEXT_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/eval/EvalContext.hpp"

namespace pdesolver::fem::eval {

	template<typename Geometry, typename Basis, Int Dim>
	class PoissonEvalElement;

	template<typename Geometry, typename Basis>
	class PoissonEvalElement<Geometry, Basis, 2> {
	public:
		static constexpr Int Dimension = 2;
		static constexpr Int NumNodes = Basis::NumNodes;
		
		// node coordinates
		const Real nodeCoords[Dimension*NumNodes];
		
		// physical coordinate
		Real x[Dimension*NumNodes];

		// time coordinate
		Real t;

		// quadrature
		Real xi[Dimension];
		Real w;

		// geometry
		Real J[Dimension*Dimension];
		Real invJ[Dimension*Dimension];
		Real detJ;

		// basis values
		Real N[NumNodes];

		// ref gradients
		Real dNdxi[NumNodes];
		Real dNdeta[NumNodes];

		// physical gradients
		Real dNdx[NumNodes];
		Real dNdy[NumNodes];

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

	}; // class PoissonEvalElement<Geometry, Basis, 2>
	
	template<typename Geometry, typename Basis>
	class PoissonEvalElement<Geometry, Basis, 3> {
	public:
		static constexpr Int Dimension = 3;
		static constexpr Int NumNodes = Basis::NumNodes;
		
		// node coordinates
		const Real nodeCoords[Dimension*NumNodes];
		
		// physical coordinate
		Real x[Dimension];

		// time coordinate
		Real t;

		// quadrature
		Real xi[Dimension];
		Real w;

		// geometry
		Real J[Dimension*Dimension];
		Real invJ[Dimension];
		Real detJ;

		// basis values
		Real N[NumNodes];

		// ref gradients
		Real dNdxi[NumNodes];
		Real dNdeta[NumNodes];

		// physical gradients
		Real dNdx[NumNodes];
		Real dNdy[NumNodes];

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

	}; // class PoissonEvalElement<Geometry, Basis, 3>


} // namespace pdesolver::fem::eval

#endif
