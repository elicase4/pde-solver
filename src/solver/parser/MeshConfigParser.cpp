#include "solver/parser/MeshConfigParser.hpp"
#include "io/YAMLReader.hpp"

pdesolver::solver::config::MeshConfig pdesolver::solver::parser::MeshConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	pdesolver::solver::config::MeshConfig cfg;

	cfg.file = YAMLReader::required<std::string>(node, "file");

	return cfg;

}
