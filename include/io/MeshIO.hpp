#ifndef PDESOLVER_MESH_IO_HPP
#define PDESOLVER_MESH_IO_HPP

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "core/Types.hpp"
#include "mesh/MeshBase.hpp"

namespace pdesolver {
	namespace io {
		
		class MeshIO {
			public:
				void writeVTK(const mesh::MeshBase& mesh, const std::string& filename);
			private:
				std::vector<Index> rowMajorToCCW(const Index* row_major_ordering, Index num_nodes);

		}; // class MeshIO
	
	} // namespace io
} // namespace pdesolver

#endif
