#ifndef PDESOLVER_IO_GMSHREADER
#define PDESOLVER_IO_GMSHREADER

#include <fstream>
#include <iostream>
#include <unordered_map>

#include "core/Types.hpp"
#include "mesh/Mesh.hpp"

namespace pdesolver {
	namespace io {
		class GmshReader {
			public:
				void read(mesh::Mesh& mesh, const std::string& filename, const std::unordered_map<int, Int>& physicalGroupMap = {});
			private:
				void readMSH4(std::istream& is, mesh::Mesh& mesh, const std::unoredered_map<int, Int>& pgMap);
				void readMSH2(std::istream& is, mesh::Mesh& mesh, const std::unoredered_map<int, Int>& pgMap);
		}
	} // namespace io
} // namespace pdesolver

#endif
