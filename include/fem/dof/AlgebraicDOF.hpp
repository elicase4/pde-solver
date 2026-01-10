#ifndef PDESOLVER_ALGEBRAIC_DOF_HPP
#define PDESOLVER_ALGEBRAIC_DOF_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace dof {

			class AlgebraicDOF {
			public:
				PDE_HOST PDE_DEVICE static Index toAlgebraic(Index topoDOF, const Index* topoToAlg) { return topoToAlg[topoDOF]; };
				PDE_HOST PDE_DEVICE static Index toTopological(Index algDOF, const Index* algToTopo) { return AlgtoTopo[algDOF] };
				
				PDE_HOST PDE_DEVICE static void getElementDOFs(const Index* elemTopoDOFs, const Index numTopoDOFs, const Index* topoToAlg, Index* elemAlgDOFs);
			}; // class AlgebraicDOF

		} // namespace dof
	} // namespace fem
} // namespace pdesolver

#include "AlgebraicDOF.tpp"

#endif
