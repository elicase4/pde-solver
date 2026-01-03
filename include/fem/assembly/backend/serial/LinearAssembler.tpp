namespace pdesolver::fem::assembly {

template<typename Element, typename Quadrature>
linalg::types::SparseMatrix<Real, Serial> LinearAssembler<Element, Quadrature, Serial>::createMatrixSystem(const mesh::MeshBase& mesh){
	
}

template<typename Element, typename Quadrature>
linalg::types::Vector<Real, Serial> LinearAssembler<Element, Quadrature, Serial>::createOperatorSystem(const mesh::MeshBase& mesh){
}

template<typename Element, typename Quadrature>
linalg::types::Vector<Real, Serial> LinearAssembler<Element, Quadrature, Serial>::createRHSVector(const mesh::MeshBase& mesh){
}

template<typename Element, typename Quadrature, typename BilinearForm>
void LinearAssembler<Element, Quadrature, Serial>::assembleMatrixSystem(const mesh::MeshBase& mesh, linalg::types::SparseMatrix<Real, Backend>& K){
}

template<typename Element, typename Quadrature, typename BilinearForm>
void LinearAssembler<Element, Quadrature, Serial>::assembleOperatorSystem(const mesh::MeshBase& mesh, linalg::types::Vector<Real, Backend>& O){
}

template<typename Element, typename Quadrature, typename LinearForm>
void LinearAssembler<Element, Quadrature, Serial>::assembleRHSVector(const mesh::MeshBase& mesh, linalg::types::Vector<Real, Backend>& F){
}

} // namespace pdesolver::fem::assembly
