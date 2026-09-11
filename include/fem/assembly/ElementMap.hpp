#ifndef PDESOLVER_FEM_ASSEMBLY_ELEMENTMAP_HPP
#define PDESOLVER_FEM_ASSEMBLY_ELEMENTMAP_HPP

#include "config/Platform.hpp"
#include "core/Types.hpp"

#include "fem/boundary/EssentialBoundaryRegistry.hpp"

#include "topology/TopologicalDOF.hpp"

namespace pdesolver {
	namespace fem {
		namespace assembly {

			enum class GatherMode { Free, Full, Constrained };

			struct NoSolution {};

			template<Index numDOFs, Index SpatialDim, GatherMode Mode, typename VectorType = NoSolution>
			PDE_HOST PDE_DEVICE void gatherElementVector(const Index* nodeIDs, Index nodesPerElement, const Real* nodeCoords, const topology::TopologicalDOF<numDOFs>& topoDOF, const fem::boundary::EssentialBoundaryRegistry* bcRegistry, const Real time, Real* Ue, const VectorType* U = nullptr);

			template<Index numDOFs, typename VectorType>
			PDE_HOST PDE_DEVICE void scatterElementVector(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Xe, VectorType& X, Real coefficient = Real(1));

			template<Index numDOFs, typename MatrixType>
			PDE_HOST PDE_DEVICE void scatterElementMatrix(const Index* nodeIDs, Index nodesPerElement, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real* Ke, MatrixType& K, Real coefficient = Real(1));

		} // namespace assembly
	} // namespace fem
} // namespace pdesolver

#include "fem/assembly/backend/cpu/ElementMap.tpp"

#endif
