namespace pdesolver::fem::assembly {

template<typename Element, typename Quadrature>
linalg::types::SparseMatrix<Real, Serial> LinearAssembler<Element, Quadrature, Serial>::createMatrixSystem(const mesh::Mesh& mesh, const Index DofsPerNode){
	
}

template<typename Element, typename Quadrature>
linalg::types::Vector<Real, Serial> LinearAssembler<Element, Quadrature, Serial>::createOperatorSystem(const mesh::Mesh& mesh, const Index DofsPerNode){
}

template<typename Element, typename Quadrature>
linalg::types::Vector<Real, Serial> LinearAssembler<Element, Quadrature, Serial>::createRHSVector(const mesh::Mesh& mesh, const Index DofsPerNode){
}

template<typename Element, typename Quadrature, typename BilinearForm>
void LinearAssembler<Element, Quadrature, Serial>::assembleMatrixSystem(const mesh::Mesh& mesh, const Index DofsPerNode, linalg::types::SparseMatrix<Real, Backend>& K){
}

template<typename Element, typename Quadrature, typename BilinearForm>
void LinearAssembler<Element, Quadrature, Serial>::assembleOperatorSystem(const mesh::Mesh& mesh, const Index DofsPerNode, linalg::types::Vector<Real, Backend>& O){
}

template<typename Element, typename Quadrature, typename LinearForm>
void LinearAssembler<Element, Quadrature, Serial>::assembleRHSVector(const mesh::Mesh& mesh, const Index DofsPerNode, linalg::types::Vector<Real, Backend>& F){
}

} // namespace pdesolver::fem::assembly
