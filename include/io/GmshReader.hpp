#ifndef PDESOLVER_IO_GMSHREADER_HPP
#define PDESOLVER_IO_GMSHREADER_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "core/Types.hpp"
#include "io/utils/Binary.hpp"
#include "mesh/Mesh.hpp"

namespace pdesolver {
	namespace io {
		class GmshReader {
			public:
				
				void read(mesh::Mesh& mesh, const std::string& filename, const std::unordered_map<int, Int>& physicalGroupMap = {});

			private:
				void readMSH4(std::istream& is, mesh::Mesh& mesh, const std::unordered_map<int, Int>& pgMap);
				void readMSH2(std::istream& is, mesh::Mesh& mesh, const std::unordered_map<int, Int>& pgMap);
		}; // class GmashReader
	} // namespace io
} // namespace pdesolver

#endif
