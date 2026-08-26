#include "application/heateq/parser/InitialConditionConfigParser.hpp"
#include "io/YAMLReader.hpp"

pdesolver::application::heateq::config::InitialConditionConfig::Type pdesolver::application::heateq::parser::InitialConditionConfigParser::parseType(const std::string& str) {

	if (str == "expression") {
		return pdesolver::application::heateq::config::InitialConditionConfig::Type::Expression;
	}

	if (str == "file") {
		return pdesolver::application::heateq::config::InitialConditionConfig::Type::File;
	}

	throw std::runtime_error("Unknown initial_condition type: " + str);

}

pdesolver::application::heateq::config::InitialConditionConfig pdesolver::application::heateq::parser::InitialConditionConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	pdesolver::application::heateq::config::InitialConditionConfig cfg;

	cfg.type = pdesolver::application::heateq::parser::InitialConditionConfigParser::parseType(YAMLReader::required<std::string>(node, "type"));

	if (cfg.type == pdesolver::application::heateq::config::InitialConditionConfig::Type::Expression) {
		cfg.expression = YAMLReader::required<std::string>(node, "expression");
	} else {
		cfg.file = YAMLReader::required<std::string>(node, "file");
	}

	return cfg;

}
