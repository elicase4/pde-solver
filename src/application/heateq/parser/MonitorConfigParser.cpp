#include "application/heateq/parser/MonitorConfigParser.hpp"

#include <stdexcept>

#include "io/YAMLReader.hpp"

pdesolver::application::heateq::config::MonitorConfig::Quantity pdesolver::application::heateq::parser::MonitorConfigParser::parseQuantity(const std::string& str) {

	if (str == "heat_flux"){
		return pdesolver::application::heateq::config::MonitorConfig::Quantity::HeatFlux;
	} else {
		throw std::runtime_error("Unknown monitor quantity: " + str);
	}

}

pdesolver::application::heateq::config::MonitorConfig pdesolver::application::heateq::parser::MonitorConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	pdesolver::application::heateq::config::MonitorConfig cfg;

	cfg.name = YAMLReader::required<std::string>(node, "name");
	cfg.quantity = pdesolver::application::heateq::parser::MonitorConfigParser::parseQuantity(YAMLReader::required<std::string>(node, "quantity"));

	const std::string reductionStr = YAMLReader::required<std::string>(node, "reduction");
	if (reductionStr == "integral") {
		cfg.reduction = pdesolver::fem::quantity::Reduction::Integral;
	} else if (reductionStr == "average") {
		cfg.reduction = pdesolver::fem::quantity::Reduction::Average;
	} else {
		throw std::runtime_error("MonitorConfigParser: monitor '" + cfg.name + "': unknown reduction '" + reductionStr + "'");
	}

	const YAML::Node& termsNode = node["terms"];
	if (!termsNode) {
		throw std::runtime_error("MonitorConfigParser: monitor '" + cfg.name + "': missing required 'terms' list");
	}

	for (const auto& termNode : termsNode) {

		pdesolver::application::heateq::config::MonitorTermConfig term;
		term.boundary = YAMLReader::required<Int>(termNode, "boundary");
		term.coefficient = YAMLReader::optional<Real>(termNode, "coefficient", Real(1.0));
		cfg.terms.push_back(term);

	}

	return cfg;

}
