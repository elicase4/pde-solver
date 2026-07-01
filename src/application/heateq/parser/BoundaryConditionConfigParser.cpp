#include "application/heateq/BoundaryConditionConfigParser.hpp"

pdesolver::application::heateq::BoundaryConditionConfig::Type pdesolver::application::heateq::BoundaryConditionConfigParser::parseBoundaryConditionType(const std::string& str) {

	if (str == "value") {
		return pdesolver::application::heateq::BoundaryConditionConfig::Type::Value;
	}

	if (str == "flux") {
		return pdesolver::application::heateq::BoundaryConditionConfig::Type::Flux;
	}

	throw std::runtime_error("Unknown boundary type: " + str);

}
