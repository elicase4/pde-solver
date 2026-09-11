#include "application/heateq/parser/ScalarMaterialPropertyConfigParser.hpp"

#include "io/YAMLReader.hpp"

pdesolver::application::heateq::config::ScalarMaterialPropertyConfig pdesolver::application::heateq::parser::ScalarMaterialPropertyConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	pdesolver::application::heateq::config::ScalarMaterialPropertyConfig cfg;

	cfg.value = YAMLReader::required<Real>(node, "value");
	cfg.unit = YAMLReader::required<std::string>(node, "unit");

	return cfg;

}
