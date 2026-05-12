#ifndef PDESOLVER_IO_UTILS_GMSH_HPP
#define PDESOLVER_IO_UTILS_GMSH_HPP

#include <array>
#include <sstream>

#include "core/Types.hpp"

namespace pdesolver {
	namespace io {
		namespace gmsh {

			struct ElemTypeInfo {
				int gmshType;
				int dim; // parametric
				Index numNodes;
				Index numFaces;
			}; // struct ElemTypeInfo

			// ---------------------------------------------------------------------------
			// Gmsh element-type table
			//
			// Current element types:
			//   Gmsh type 3  → 4-node quad  (2-D, nodesPerElement = 4, facesPerElement = 4)
			//   Gmsh type 5  → 8-node hex   (3-D, nodesPerElement = 8, facesPerElement = 6)
			//   Gmsh type 1  → 2-node line  (boundary / ignored as volumetric)
			//   Gmsh type 15 → 1-node point (ignored as volumetric)
			//
			// ---------------------------------------------------------------------------
			static constexpr std::array<ElemTypeInfo, 4> ELEM_TABLE = {{
				{3, 2, 4, 4}, // quad4
				{5, 3, 8, 6}, // hex8
				{1, 1, 2, 0}, // line2
				{15, 0, 1, 0}, // point
			}};

			static const ElemTypeInfo* lookupElemType(int gmshType) {
				for (const auto& e : ELEM_TABLE){
					if (e.gmshType == gmshType) return &e;
				}
				return nullptr;
			}

			// skip whitespace
			static bool skipWS(std::istream& is){
				while ((is.peek() == ' ') || (is.peek() == '\t') || (is.peek() == '\n') || (is.peek() == '\r')){
					is.get();
				}
				return !is.eof();
			}

			// read an entire line, skipping any blank lines
			static std::string nextLine(std::istream& is){
				
				std::string line;

				while (std::getline(is, line)) {
					// skip trailing cursor return
					if (!line.empty() && line.back() == '\r') line.pop_back();
					if (!line.empty()) return line;
				}

				return {};
			}

		} // namespace gmsh
	} // namespace io
} // namespace pdesolver

#endif
