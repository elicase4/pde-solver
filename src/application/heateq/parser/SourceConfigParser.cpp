#include "application/heateq/parser/SourceConfigParser.hpp"

#include <stdexcept>

#include "io/YAMLReader.hpp"
#include "solver/parser/NodalFieldReadConfigParser.hpp"

pdesolver::application::heateq::config::SourceConfig::Type pdesolver::application::heateq::parser::SourceConfigParser::parseSourceType(const std::string& str) {

	if (str == "volumetric") {
		return pdesolver::application::heateq::config::SourceConfig::Type::VolumetricHeatSource;
	}

	throw std::runtime_error("Unknown source type: " + str);

}

pdesolver::application::heateq::config::SourceConfig pdesolver::application::heateq::parser::SourceConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	pdesolver::application::heateq::config::SourceConfig cfg;

	cfg.type = pdesolver::application::heateq::parser::SourceConfigParser::parseSourceType(YAMLReader::required<std::string>(node, "type"));

	const YAML::Node& readNode = node["read"];
	if (!readNode) {
		throw std::runtime_error("SourceConfigReader: missing required 'read' section in 'physics.models.source'");
	}

	cfg.read = pdesolver::solver::parser::NodalFieldReadConfigParser::parse(readNode);

	return cfg;

}
