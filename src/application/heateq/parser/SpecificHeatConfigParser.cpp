#include "application/heateq/parser/SpecificHeatConfigParser.hpp"
#include "io/YAMLReader.hpp"

pdesolver::application::heateq::config::SpecificHeatConfig::Type pdesolver::application::heateq::parser::SpecificHeatConfigParser::parseSpecificHeatType(const std::string& str) {

	if (str == "constant") {
		return pdesolver::application::heateq::config::SpecificHeatConfig::Type::Constant;
	}

	if (str == "temperature_dependent") {
		return pdesolver::application::heateq::config::SpecificHeatConfig::Type::TemperatureDependent;
	}

	throw std::runtime_error("Unknown specific heat type: " + str);

}

pdesolver::application::heateq::config::SpecificHeatConfig pdesolver::application::heateq::parser::SpecificHeatConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;
	using Type = pdesolver::application::heateq::config::SpecificHeatConfig::Type;

	pdesolver::application::heateq::config::SpecificHeatConfig cfg;

	cfg.type = node["type"] ? pdesolver::application::heateq::parser::SpecificHeatConfigParser::parseSpecificHeatType(YAMLReader::required<std::string>(node, "type")) : Type::Constant;

	if (cfg.type == Type::Constant) {
		cfg.value = YAMLReader::required<Real>(node, "value");
	} else if (cfg.type == Type::TemperatureDependent) {
		cfg.valueExpression = YAMLReader::required<std::string>(node, "value_expression");
		cfg.gradientExpression = YAMLReader::required<std::string>(node, "gradient_expression");
	}

	cfg.unit = YAMLReader::required<std::string>(node, "unit");

	return cfg;

}
