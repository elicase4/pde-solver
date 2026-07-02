#include "solver/config/MeshConfigParser.hpp"
#include "io/YAML/Reader.hpp"

pdesolver::solver::config::MeshConfig pdesolver::solver::parser::MeshConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	MeshConfig cfg;

	cfg.file = YAMLReader::required<Index>(node, "file");

	return cfg;

}
