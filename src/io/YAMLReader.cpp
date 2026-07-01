#include "io/YAMLReader.hpp"

YAML::Node pdesolver::io::YAMLReader::loadFile(const std::string& filename) {
	return YAML::LoadFile(filename);
}
