#include "application/heateq/SourceConfigParser.hpp"
#include "io/YAMLReader.hpp"

pdesolver::application::heateq::config::SourceConfig::Type pdesolver::application::heateq::parser::SourceConfigParser::parseSourceType(const std::string& str) {

	if (str == "volumetric_heat_source") {
		return pdesolver::application::heateq::ConductivityConfig::Type::VolumetricHeatSource;
	}

	throw std::runtime_error("Unknown source type: " + str);

}

pdesolver::application::heateq::config::SourceConfig pdesolver::application::heateq::parser::SourceConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	SourceConfig cfg;

	cfg.source.expression = YAMLReader::required<std::string>(src, "expression");
	
	return cfg;

}
