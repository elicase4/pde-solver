#include "io/GmshReader.hpp"

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

void pdesolver::io::GmshReader::read(mesh::Mesh& mesh, const std::string& filename, const std::unordered_map<int, Int>& physicalGroupMap){

}

void pdesolver::io::GmshReader::readMSH4(std::istream& is, mesh::Mesh& mesh, const std::unordered_map<int, Int>& pgMap){

}

void pdesolver::io::GmshReader::readMSH2(std::istream& is, mesh::Mesh& mesh, const std::unordered_map<int, Int>& pgMap){

}
