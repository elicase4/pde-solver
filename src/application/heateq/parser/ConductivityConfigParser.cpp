#include "application/heateq/ConductivityConfigParser.hpp"

pdesolver::application::heateq::ConductivityConfig::Type pdesolver::application::heateq::ConductivityConfigParser::parseConductivityType(const std::string& str) {

	if (str == "constant") {
		return pdesolver::application::heateq::ConductivityConfig::Type::Constant;
	}

	if (str == "anisotropic") {
		return pdesolver::application::heateq::ConductivityConfig::Type::Anisotropic;
	}

	throw std::runtime_error("Unknown conductivity type: " + str);

}
