#ifndef HEATEQUATION_MASSFORM_HPP
#define HEATEQUATION_MASSFORM_HPP

#include "fem/form/BilinearForm.hpp"

namespace pdesolver::equations::heateq {

	template<typename QuadraturePointVolume>
	struct MassForm {

		PDE_HOST PDE_DEVICE static void computeElementLevelMatrix(const QuadraturePointVolume& qp, const Real*, Real* Me){
			
			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index b = 0; b < qp.nodesPerElement(); ++b){
					Me[a * qp.nodesPerElement() + b] += qp.N[a] * qp.N[b] * qp.measure * qp.w;
				}
			}
		}
		
		PDE_HOST PDE_DEVICE static void computeElementLevelVector(const QuadraturePointVolume& qp, const Real* Ue, Real* Oe){
			
			for (Index a = 0; a < qp.nodesPerElement(); ++a){
				for (Index b = 0; b < qp.nodesPerElement(); ++b){
					
					Oe[a] += qp.N[a] * qp.N[b] * Ue[b] * qp.measure * qp.w;
				}
			}
		}

	}; // struct MassForm

} // namespace pdesolver::equations::heateq

#endif
