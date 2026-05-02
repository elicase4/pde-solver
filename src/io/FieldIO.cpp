#include "io/FieldIO.hpp"

void pdesolver::io::FieldIO::writeVTK(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const fem::boundary::BoundaryRegistry& bcRegistry, Real time, const Real* algField, const std::vector<std::string>& dofNames, const std::string& filename, VTKWriter::Format fmt = VTKWriter::Format::ASCII){

}

std::vector<Real> pdesolver::io::FieldIO::reconstructNodalField(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const fem::boundary::BoundaryRegistry& bcRegistry, Real time, const Real* algField) const {

}
