#ifndef PDESOLVER_IO_GMSHREADER_HPP
#define PDESOLVER_IO_GMSHREADER_HPP

#include <cstdint>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "core/Types.hpp"
#include "io/utils/Binary.hpp"
#include "io/utils/Gmsh.hpp"
#include "mesh/exchange/gmsh/IntermediateMesh.hpp"

namespace pdesolver {
	namespace io {
		class GmshReader {
			public:
				
				static void read(mesh::exchange::gmsh::IntermediateMesh& mesh, const std::string& filename);

			private:
				
				enum class Format {
					ASCII,
					Binary
				}; // enum class Format

				struct VersionInfo {
					double version;
					Format format;
				};
				
				// high level version readers
				static VersionInfo readMeshFormat(std::istream& is);

				static void readMSH4(std::istream& is, mesh::exchange::gmsh::IntermediateMesh& mesh, VersionInfo::Format fmt);
				
				static void readMSH2(std::istream& is, mesh::exchange::gmsh::IntermediateMesh& mesh, VersionInfo::Format fmt);

				// internal helpers
				static void skipToSection(std::istream& is, const std::string& tag);

				static void deduceDimensions(mesh::exchange::gmsh::IntermediateMesh& mesh);

				static void readPhysicalNames(std::istream& is, std::unordered_map<Int, std::string>& names);

				static std::unordered_map<Int, Int> readEntities(std::istream& is, VersionInfo::Format fmt);

				static void readNodes(std::istream& is, mesh::exchange::gmsh::IntermediateMesh& mesh, std::unordered_map<Index, Index>& tagToIdx, VersionInfo::Format fmt);
				
				static void readElements(std::istream& is, mesh::exchange::gmsh::IntermediateMesh& mesh, const std::unordered_map<Index, Index>& tagToIdx, const std::unordered_map<Int, Int>& entityPhys, VersionInfo::Format fmt);

		}; // class GmshReader
	} // namespace io
} // namespace pdesolver

#endif
