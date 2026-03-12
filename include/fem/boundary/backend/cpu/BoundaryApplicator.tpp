namespace pdesolver::fem::boundary {

void BoundaryApplicator<linalg::types::backend::CPU>::applyEssentialBC(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& registry, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

}

void BoundaryApplicator<linalg::types::backend::CPU>::applyNaturalBC(const mesh::Mesh& mesh, const topology::TopologicalDOF& topoDOF, const BoundaryRegistry& registry, linalg::types::Vector<Real, linalg::types::backend::CPU>& F){

}

} // namespace pdesolver::fem::boundary
