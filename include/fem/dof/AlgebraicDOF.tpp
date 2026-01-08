namespace pdesolver::fem::dof {

PDE_HOST PDE_DEVICE void AlgebraicDOF::getElementDOFs(const Index* elemTopoDOFs, Index numTopoDOFs, const Index* topoToAlg, Index* elemAlgDOFs){
	
	Index count = 0;
	
	for (Index i = 0; i < numTopoDOFs; ++i) {
		Index algDOF = topoToAlg[elemTopoDOFs[i]];
		if (algDOF != -1){
			elemAlgDOFs[count++] = algDOF;
		}
	}

	return count;

}

} // namespace pdesolver::fem::dof
