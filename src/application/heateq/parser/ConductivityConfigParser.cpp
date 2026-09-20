#include "application/heateq/parser/ConductivityConfigParser.hpp"
#include "io/YAMLReader.hpp"

pdesolver::application::heateq::config::ConductivityConfig::Type pdesolver::application::heateq::parser::ConductivityConfigParser::parseConductivityType(const std::string& str) {

	if (str == "constant") {
		return pdesolver::application::heateq::config::ConductivityConfig::Type::Constant;
	}

	if (str == "anisotropic") {
		return pdesolver::application::heateq::config::ConductivityConfig::Type::Anisotropic;
	}

	if (str == "temperature_dependent_isotropic") {
		return pdesolver::application::heateq::config::ConductivityConfig::Type::TemperatureDependentIsotropic;
	}

	if (str == "temperature_dependent_anisotropic") {
		return pdesolver::application::heateq::config::ConductivityConfig::Type::TemperatureDependentAnisotropic;
	}

	throw std::runtime_error("Unknown conductivity type: " + str);

}

pdesolver::application::heateq::config::ConductivityConfig pdesolver::application::heateq::parser::ConductivityConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	pdesolver::application::heateq::config::ConductivityConfig cfg;

	cfg.type = pdesolver::application::heateq::parser::ConductivityConfigParser::parseConductivityType(YAMLReader::required<std::string>(node, "type"));

	using Type = pdesolver::application::heateq::config::ConductivityConfig::Type;

	if (cfg.type == Type::Constant) {
		cfg.value = YAMLReader::required<Real>(node, "value");
	} else if (cfg.type == Type::Anisotropic) {
		cfg.tensor = YAMLReader::required<std::vector<std::vector<Real>>>(node, "tensor");
	} else if (cfg.type == Type::TemperatureDependentIsotropic) {
		cfg.valueExpression = YAMLReader::required<std::string>(node, "value_expression");
		cfg.gradientExpression = YAMLReader::required<std::string>(node, "gradient_expression");
	} else if (cfg.type == Type::TemperatureDependentAnisotropic) {
		cfg.tensor = YAMLReader::required<std::vector<std::vector<Real>>>(node, "tensor");
		cfg.valueExpression = YAMLReader::required<std::string>(node, "value_expression");
		cfg.gradientExpression = YAMLReader::required<std::string>(node, "gradient_expression");
	}

	cfg.unit = YAMLReader::required<std::string>(node, "unit");

	return cfg;

}
