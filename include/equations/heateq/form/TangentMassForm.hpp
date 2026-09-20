#ifndef HEATEQUATION_TANGENTMASSFORM_HPP
#define HEATEQUATION_TANGENTMASSFORM_HPP

#include "fem/form/BilinearForm.hpp"

namespace pdesolver::equations::heateq {

	template<typename QuadraturePointVolume>
	struct TangentMassForm {

		PDE_HOST PDE_DEVICE static void computeElementLevelMatrix(const QuadraturePointVolume& qp, Real* Me){

			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index b = 0; b < qp.nodesPerElement(); ++b){
					Me[a * qp.nodesPerElement() + b] += qp.rho * qp.dcpdT * qp.T.rate * qp.N[a] * qp.N[b] * qp.measure * qp.w;
				}
			}
		}

		PDE_HOST PDE_DEVICE static void computeElementLevelVector(const QuadraturePointVolume& qp, const Real* Ue, Real* Oe){

			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index b = 0; b < qp.nodesPerElement(); ++b){
					Oe[a] += qp.rho * qp.dcpdT * qp.T.rate * qp.N[a] * qp.N[b] * Ue[b] * qp.measure * qp.w;
				}
			}
		}

	}; // struct TangentMassForm

} // namespace pdesolver::equations::heateq

#endif
