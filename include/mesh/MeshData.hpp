#ifndef MESH_DATA_HPP
#define MESH_DATA_HPP

#include "core/Types.hpp"
#include <vector>

namespace pdesolver {
	
	struct MeshData {
		
		// dimension
		Index parametricDim;
		Index spatialDim;
		
		// order
		std::vector<Index> basisOrder;
		
		// count
		Index numNodes;
		Index numElements;
		Index numBoundaries;
		Index nodesPerElement;
		Index nodesPerBoundaryEntity;

		// mesh data
		std::vector<Real> xyz;
		std::vector<Index> ien;
		std::vector<Index> rng;
		std::vector<Real> C;
		
		// helper methods
		bool isIGA() const { return !C.empty(); }
		
		Real* getNodeCoord(Index nodeID) {
			return &xyz[nodeID * spatialDim];
		}
		
		Index* getElementNodes(Index elemID) {
			return &ien[elemID *nodesPerElement];
		}

		Int getBoundaryTag(Index nodeID) const {
			return rng[nodeID];
		}

		void setBoundaryTag(Index nodeID, Int tag) {
			rng[nodeID] = tag;
		}
		
		bool isOnBoundary(Index nodeID) const {
			return rng[nodeID] >= 0;
		}

		bool isOnBoundaryTag(Index NodeID, Int tag){
			return rng[nodeID] == tag;
		}
		
		Index extractionMatrixSize() const;

		Real* getExtractionOperator(Index elemID);
		
		void clear() {
			parametricDim = 0; spatialDim = 0; basisOrder = 0; numNodes = 0; numElements = 0; nodesPerElement = 0; nodesPerBoundaryEntity = 0;
			xyz.clear(); ien.clear(); rng.clear(); C.clear();
		}

		bool isValid() const {
			
			if (xyz.size() != numNodes * spatialDim) {
				return false;
			}
			if (ien.size() != numElements * nodesPerElement) {
				return false;
			}
			if (rng.size() != numNodes) {
				return false;
			}
			
			if (parametricDim == 0 || spatialDim == 0) {
				return false;
			}
			
			for (Index p : basisOrder){
				if (p == 0){
					return false;
				}
			}
			
			if (nodesPerElement != computeNodesPerElement(basisOrder)) {
				return false;
			}
			
			return true;
		}
		
		static Index computeNodesPerElement(const std::vector<Index>& basisOrder) {
			Index npe = 1;
			for (Index p : basisOrder){
				npe *= (p + 1);
			}
			return npe;
		}

	}

}

#endif
