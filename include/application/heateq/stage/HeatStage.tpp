namespace pdesolver::application::heateq::stage {

	template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
	HeatStage<BackendType, BasisType, QuadratureVolumeType, QuadratureBoundaryType>::HeatStage(const HeatConfig& config) : config_(config) {}

	template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
	void HeatStage<BackendType, BasisType, QuadratureVolumeType, QuadratureBoundaryType>::initialize() {

	}

	template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
	void HeatStage<BackendType, BasisType, QuadratureVolumeType, QuadratureBoundaryType>::assemble() {

	}

	template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
	bool HeatStage<BackendType, BasisType, QuadratureVolumeType, QuadratureBoundaryType>::solve() {
		
		return true;
	
	}
	
	template<typename BackendType, typename BasisType, typename QuadratureVolumeType, typename QuadratureBoundaryType>
	void HeatStage<BackendType, BasisType, QuadratureVolumeType, QuadratureBoundaryType>::finalize() {

	}

} // namespace pdesolver::application::heateq::stage
