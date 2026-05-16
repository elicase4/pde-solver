#ifndef PDESOLVER_IO_GMSHREADER_HPP
#define PDESOLVER_IO_GMSHREADER_HPP

#include <istream>
#include <string>

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

				static VersionInfo readMeshFormat(std::istream& is);

				static void readMSH2ASCII(std::istream& is, mesh::exchange::gmsh::IntermediateMesh& mesh);
				
				static void readMSH2Binary(std::istream& is, mesh::exchange::gmsh::IntermediateMesh& mesh);

				static void readMSH4ASCII(std::istream& is, mesh::exchange::gmsh::IntermediateMesh& mesh);
				
				static void readMSH4Binary(std::istream& is, mesh::exchange::gmsh::IntermediateMesh& mesh);

		}; // class GmshReader
	} // namespace io
} // namespace pdesolver

#endif
