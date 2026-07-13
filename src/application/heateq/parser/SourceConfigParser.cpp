#include "application/heateq/parser/SourceConfigParser.hpp"
#include "io/YAMLReader.hpp"

pdesolver::application::heateq::config::SourceConfig::Type pdesolver::application::heateq::parser::SourceConfigParser::parseSourceType(const std::string& str) {

	if (str == "volumetric_heat_source") {
		return pdesolver::application::heateq::config::SourceConfig::Type::VolumetricHeatSource;
	}

	throw std::runtime_error("Unknown source type: " + str);

}

pdesolver::application::heateq::config::SourceConfig pdesolver::application::heateq::parser::SourceConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	pdesolver::application::heateq::config::SourceConfig cfg;

	cfg.expression = YAMLReader::required<std::string>(node, "expression");
	
	return cfg;

}
