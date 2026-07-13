#include "application/heateq/parser/ConductivityConfigParser.hpp"
#include "io/YAMLReader.hpp"

pdesolver::application::heateq::config::ConductivityConfig::Type pdesolver::application::heateq::parser::ConductivityConfigParser::parseConductivityType(const std::string& str) {

	if (str == "constant") {
		return pdesolver::application::heateq::config::ConductivityConfig::Type::Constant;
	}

	if (str == "anisotropic") {
		return pdesolver::application::heateq::config::ConductivityConfig::Type::Anisotropic;
	}

	throw std::runtime_error("Unknown conductivity type: " + str);

}

pdesolver::application::heateq::config::ConductivityConfig pdesolver::application::heateq::parser::ConductivityConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	pdesolver::application::heateq::config::ConductivityConfig cfg;

	cfg.type = pdesolver::application::heateq::parser::ConductivityConfigParser::parseConductivityType(YAMLReader::required<std::string>(node, "type"));

	if (cfg.type == pdesolver::application::heateq::config::ConductivityConfig::Type::Constant) {
		cfg.value = YAMLReader::required<Real>(node, "value");
	} else if (cfg.type == pdesolver::application::heateq::config::ConductivityConfig::Type::Anisotropic) {
		cfg.tensor = YAMLReader::required<std::vector<std::vector<Real>>>(node, "tensor");
	}

	return cfg;

}
