namespace pdesolver::fem::assembly {

linalg::types::SparseMatrix<Real, Backend> createMatrixSystem(const mesh::MeshBase& mesh, const topology::TopologicalDOF& topoDOF){
}

linalg::types::Vector<Real, Backend> createOperatorSystem(const mesh::MeshBase& mesh, const topology::TopologicalDOF& topoDOF){
}

linalg::types::Vector<Real, Backend> createRHSVector(const mesh::MeshBase& mesh, const topology::TopologicalDOF& topoDOF){
}

template<typename BilinearForm>
void assembleMatrixSystem(const mesh::MeshBase& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::SparseMatrix<Real, Backend>& K){
}

template<typename BilinearForm>
void assembleOperatorSystem(const mesh::MeshBase& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::Vector<Real, Backend>& O){
}

template<typename LinearForm>
void assembleRHSVector(const mesh::MeshBase& mesh, const topology::TopologicalDOF& topoDOF, linalg::types::Vector<Real, Backend>& F){
}

} // namespace pdesolver::fem::assembly
