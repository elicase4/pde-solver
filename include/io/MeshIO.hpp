#ifndef PDESOLVER_MESH_IO_HPP
#define PDESOLVER_MESH_IO_HPP

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "core/Types.hpp"
#include "mesh/Mesh.hpp"
#include "VTKWriter.hpp"

/* Binary Mesh Format Description
[HEADER]
  magic:         uint32  = 0x504D5348  ("PMSH")
  version:       uint32  = 1
  parametricDim: uint32
  spatialDim:    uint32
  numNodes:      uint64
  numElements:   uint64
  nodesPerElem:  uint32
  facesPerElem:  uint32
  numBasisOrders:uint32
  basisOrders:   uint32[numBasisOrders]
  hasExtractionOps: uint8  (0 = FEM, 1 = IGA)
  extractionOpSize: uint64 (0 if FEM)

[DATA]
  xyz:  Real[numNodes * spatialDim]       (floats/doubles, little-endian)
  ien:  Index[numElements * nodesPerElem] (uint64, little-endian)
  rng:  Int[numElements * facesPerElem]   (int32, little-endian)
  C:    Real[extractionOpSize]            (optional, IGA only)
*/

namespace pdesolver {
	namespace io {
		
		class MeshIO {
		public:
			void writeBinary(const mesh::Mesh& mesh, const std::string& filename);
			void readBinary(mesh::Mesh& mesh, const std::string& filename);
		}; // class MeshIO
	
	} // namespace io
} // namespace pdesolver

#endif
