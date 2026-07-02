#include "application/heateq/ConductivityConfigParser.hpp"
#include "io/YAMLReader.hpp"

pdesolver::application::heateq::config::ConductivityConfig::Type pdesolver::application::heateq::parser::ConductivityConfigParser::parseConductivityType(const std::string& str) {

	if (str == "constant") {
		return pdesolver::application::heateq::ConductivityConfig::Type::Constant;
	}

	if (str == "anisotropic") {
		return pdesolver::application::heateq::ConductivityConfig::Type::Anisotropic;
	}

	throw std::runtime_error("Unknown conductivity type: " + str);

}

pdesolver::application::heateq::config::ConductivityConfig pdesolver::application::heateq::parser::ConductivityConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	ConductivityConfig cfg;

	cfg.conductivity.type = ConductivityConfigParser::parseConductivityType(YAMLReader::required<std::string>(cond, "type"));

	if (cfg.conductivity.type == ConductivityConfig::Type::Constant) {
		cfg.conductivity.value = YAMLReader::required<Real>(cond, "value");
	} else if (cfg.conductivity.type == ConductivityConfig::Type::Anisotropic) {
		cfg.conductivity.tensor = YAMLReader::required<std::vector<std::vector<Real>>>(cond, "tensor");
	}

	return cfg;

}
