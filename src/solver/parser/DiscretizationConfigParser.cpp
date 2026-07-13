#include "solver/parser/DiscretizationConfigParser.hpp"
#include "io/YAMLReader.hpp"

pdesolver::solver::config::DiscretizationConfig pdesolver::solver::parser::DiscretizationConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	config::DiscretizationConfig cfg;

	cfg.quadraturePointXi = YAMLReader::optional<Index>(node, "quadrature_points_xi", 2);
	cfg.quadraturePointEta = YAMLReader::optional<Index>(node, "quadrature_points_eta", 2);
	cfg.quadraturePointZeta = YAMLReader::optional<Index>(node, "quadrature_points_zeta", 2);
	cfg.blockDOFOrdering = YAMLReader::optional<bool>(node, "block_dof_ordering", true);

	return cfg;

}
