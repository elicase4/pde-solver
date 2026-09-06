#include "solver/parser/NodalFieldWriteConfigParser.hpp"

#include "io/YAMLReader.hpp"

pdesolver::solver::config::NodalFieldWriteConfig pdesolver::solver::parser::NodalFieldWriteConfigParser::parse(const YAML::Node& node) {

	pdesolver::solver::config::NodalFieldWriteConfig cfg;

	const YAML::Node& writeNode = node["write"];
	if (writeNode) {
		cfg.file = pdesolver::io::YAMLReader::required<std::string>(writeNode, "file");
	}

	return cfg;

}
