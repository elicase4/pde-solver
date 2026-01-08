#ifndef PDESOLVER_ALGEBRAIC_DOF_HPP
#define PDESOLVER_ALGEBRAIC_DOF_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace pdesolver {
	namespace fem {
		namespace dof {

			class AlgebraicDOF {
			public:
				PDE_HOST PDE_DEVICE static Index toAlgebraic(Index topoDOF, const Index* topoToAlg);
				PDE_HOST PDE_DEVICE static Index toTopological(Index algDOF, const Index* algToTopo);
				PDE_HOST PDE_DEVICE static void getElementDOFs(const Index* elemTopoDOFs, Index numTopoDOFs, const Index* topoToAlg, Index* elemAlgDOFs);
			}; // class AlgebraicDOF

		} // namespace dof
	} // namespace fem
} // namespace pdesolver

#include "ALgebraicDOF.tpp"

#endif
