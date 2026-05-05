#include "io/GmshReader.hpp"


void pdesolver::io::GmshReader::read(mesh::Mesh& mesh, const std::string& filename, const std::unordered_map<int, Int>& physicalGroupMap){

	std::ifstream ifs(filename);
	if (!ifs.is_open()){
		throw std::runtime_error("GmshReader: cannot open file " + filename);
	}

	// find MeshFormat to detect version
	std::string line;
	while (std::getline(ifs, line)){
		if (line.find("$MeshFormat") != std::string::npos) break;
	}

	if (ifs.eof()){
		throw std::runtime_error("GmshReader: missing $MeshFormat section");
	}

	double version = 0.0;
	int filetype = 0;
	int datasize = 0;
	ifs >> version >> filetype >> datasize;

	if (filetype != 0){
		throw std::runtime_error("GmshReader: only ASCII Gmsh files supported");
	}

	// rewind input stream
	ifs.seekg(0);
	mesh.clear();

	// dispatch to each version reader
	if (version >= 4.0){
		readMSH4(ifs, mesh, physicalGroupMap);
	} else {
		readMSH2(ifs, mesh, physicalGroupMap);
	}

	if (!mesh.isValid()){
		throw std::runtime_error("GmshReader: mesh read from " + filename + " failed validation");
	}

	std::cout << "GmshReader: read " << mesh.data.numNodes << " nodes, " << mesh.data.numElements << " elements from " << filename << "\n";

}

void pdesolver::io::GmshReader::readMSH2(std::istream& is, mesh::Mesh& mesh, const std::unordered_map<int, Int>& pgMap){

	// Nodes Section
	{
		
		// find nodes section
		std::string line;
		while (std::getline(is, line)){
			if (line.find("$Nodes") != std::string::npos) break;
		}
		if (is.eof()){
			throw std::runtime_error("GmshReader MSH2: missing $Nodes section");
		}

		// parameters
		Index numNodes = 0;
		is >> numNodes;
		mesh.data.numNodes = numNodes;
		mesh.data.spatialDim = 3; // Gmsh always has 3D, trimming is done to go to 2D
		mesh.data.xyz.resize(numNodes*3);

		// get node coordinates
		for (Index n = 0; n < numNodes; ++n) {
			int nodeTag;
			double x, y, z;
			is >> nodeTag >> x >> y >> z;
			mesh.data.xyz[n*3+0] = static_cast<Real>(x);
			mesh.data.xyz[n*3+1] = static_cast<Real>(y);
			mesh.data.xyz[n*3+2] = static_cast<Real>(z);
		}

	}

	// Elements and boundary entities section
	{

		// find elements section
		std::string line;
		while (std::getline(is, line)){
			if (line.find("$Elements") != std::string::npos) break;
		}
		if (is.eof()){
			throw std::runtime_error("GmshReader MSH2: missing $Elements section");
		}

		// parameters
		Index totalElems = 0;
		if >> totalElems;
		std::getline(is, line); // get remaining count line

		// temporary storage to separate volume elemnts from boundary edges and faces
		struct RawElem {
			int gmshType;
			int physTag;
			std::vector<int> nodeIds;
		};

		std::vector<rawElem> volumeElems;
		std::vector<rawElem> boundaryElems;

		int detectedDim = 0;

		for (Index i = 0; i < totalElems; ++i){

			int elemTag, gmshType, numTags;
			is >> elemTag >> gmshType >> numTags;

           int physTag = 0;
            for (int t = 0; t < numTags; ++t) {
                int tag; is >> tag;
                if (t == 0) physTag = tag; // first tag = physical group
            }
 
            const ElemTypeInfo* info = lookupElemType(gmshType);
            if (!info) {
                // Skip unknown element types: consume the rest of the line
                std::getline(is, line);
                continue;
            }
 
            std::vector<int> nodes(info->numNodes);
            for (Index k = 0; k < info->numNodes; ++k)
                is >> nodes[k];
 
            if (info->dim >= 2) {
                detectedDim = std::max(detectedDim, info->dim);
                volumeElems.push_back({ gmshType, physTag, std::move(nodes) });
            } else if (info->dim == 1) {
                boundaryElems.push_back({ gmshType, physTag, std::move(nodes) });
            }
            // dim == 0 (points) are silently ignored
        }
 
        if (volumeElems.empty())
            throw std::runtime_error("GmshReader MSH2: no volume elements found");
 
        // Validate homogeneity
        for (const auto& e : volumeElems) {
            if (e.gmshType != volumeElems[0].gmshType)
                throw std::runtime_error(
                    "GmshReader MSH2: mixed element types — only homogeneous meshes are supported");
        }
 
        const ElemTypeInfo* typeInfo = lookupElemType(volumeElems[0].gmshType);
 
        d.parametricDim    = static_cast<Index>(detectedDim);
        d.spatialDim       = static_cast<Index>(detectedDim); // trim to dim
        d.numElements      = volumeElems.size();
        d.nodesPerElement  = typeInfo->numNodes;
        d.facesPerElement  = typeInfo->numFaces;
        d.basisOrder       = std::vector<Index>(detectedDim, 1); // linear
 
        // Trim xyz to spatialDim (remove z-column for 2-D meshes)
        if (d.spatialDim == 2) {
            std::vector<Real> xyz2(d.numNodes * 2);
            for (Index n = 0; n < d.numNodes; ++n) {
                xyz2[n * 2 + 0] = d.xyz[n * 3 + 0];
                xyz2[n * 2 + 1] = d.xyz[n * 3 + 1];
            }
            d.xyz = std::move(xyz2);
        }
 
        // Build a node-tag → 0-based index map (MSH2 tags are 1-based)
        // Gmsh guarantees sequential tags starting at 1 for simple meshes,
        // but we handle the general case.
        std::unordered_map<int, Index> nodeTagToIdx;
        nodeTagToIdx.reserve(d.numNodes);
        {
            // Re-read section to get the original tags — cheaper than a full
            // second pass: we stored coords already, so just rebuild the map
            // from node order (MSH2 tags are listed sequentially in $Nodes).
            // We do this by re-reading the $Nodes block.
            is.clear();
            is.seekg(0);
            std::string tmp;
            while (std::getline(is, tmp))
                if (tmp.find("$Nodes") != std::string::npos) break;
            Index nn; is >> nn;
            for (Index n = 0; n < nn; ++n) {
                int tag; double x, y, z;
                is >> tag >> x >> y >> z;
                nodeTagToIdx[tag] = n;
            }
        }
 
        // IEN
        d.ien.resize(d.numElements * d.nodesPerElement);
        for (Index e = 0; e < d.numElements; ++e) {
            for (Index k = 0; k < d.nodesPerElement; ++k) {
                int gmshTag = volumeElems[e].nodeIds[k];
                auto it = nodeTagToIdx.find(gmshTag);
                if (it == nodeTagToIdx.end())
                    throw std::runtime_error(
                        "GmshReader MSH2: unknown node tag " +
                        std::to_string(gmshTag));
                d.ien[e * d.nodesPerElement + k] = it->second;
            }
        }
 
        // RNG (boundary tags) — initialise to -1 (interior)
        d.rng.assign(d.numElements * d.facesPerElement, -1);
 
        // For each boundary edge/face, find which volume element face it belongs to
        // and record the physical tag.
        //
        // Strategy: build a map from sorted node-set → (elemIdx, faceIdx)
        // for each face of every volume element, then look up boundary elements.
        //
        // Face connectivity per element type (local node indices, 0-based):
        //   quad4 faces (edges): {0,1}, {1,2}, {2,3}, {3,0}  (CCW in Gmsh)
        //   hex8  faces:         {0,1,2,3},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}
        //
        // Note: these are Gmsh MSH2 local node orderings (NOT pdesolver CCW),
        // but since we use sorted node sets for matching, ordering doesn't matter.
 
        using NodeSet = std::vector<Index>;
        auto makeKey = [](NodeSet s) -> NodeSet {
            std::sort(s.begin(), s.end());
            return s;
        };
 
        // face definitions: [faceIdx] = list of local node indices
        std::vector<std::vector<Index>> faceDefs;
        if (d.nodesPerElement == 4) {
            faceDefs = { {0,1}, {1,2}, {2,3}, {3,0} };
        } else if (d.nodesPerElement == 8) {
            faceDefs = { {0,1,2,3}, {4,5,6,7}, {0,1,5,4},
                         {1,2,6,5}, {2,3,7,6}, {3,0,4,7} };
        }
 
        // Map sorted-node-set → (elemIdx, faceIdx)
        struct FaceRef { Index elem; Index face; };
        std::unordered_map<std::string, FaceRef> faceMap;
        faceMap.reserve(d.numElements * d.facesPerElement);
 
        auto nodeSetToKey = [](const NodeSet& s) -> std::string {
            std::string k;
            k.reserve(s.size() * sizeof(Index));
            for (Index v : s)
                k.append(reinterpret_cast<const char*>(&v), sizeof(Index));
            return k;
        };
 
        for (Index e = 0; e < d.numElements; ++e) {
            for (Index f = 0; f < d.facesPerElement; ++f) {
                NodeSet ns;
                ns.reserve(faceDefs[f].size());
                for (Index li : faceDefs[f])
                    ns.push_back(d.ien[e * d.nodesPerElement + li]);
                auto key = nodeSetToKey(makeKey(ns));
                faceMap[key] = { e, f };
            }
        }
 
        // Walk boundary elements and stamp RNG
        for (const auto& be : boundaryElems) {
            NodeSet ns;
            ns.reserve(be.nodeIds.size());
            for (int gmshTag : be.nodeIds) {
                auto it = nodeTagToIdx.find(gmshTag);
                if (it == nodeTagToIdx.end())
                    throw std::runtime_error(
                        "GmshReader MSH2: unknown boundary node tag " +
                        std::to_string(gmshTag));
                ns.push_back(it->second);
            }
            auto key = nodeSetToKey(makeKey(ns));
            auto fit = faceMap.find(key);
            if (fit == faceMap.end()) continue; // not on a known face
 
            Int solverTag = be.physTag;
            if (!pgMap.empty()) {
                auto pit = pgMap.find(be.physTag);
                if (pit != pgMap.end()) solverTag = pit->second;
            }
            d.rng[fit->second.elem * d.facesPerElement + fit->second.face] = solverTag;
        }
    }
}
 
// ---------------------------------------------------------------------------
// MSH4 reader
// ---------------------------------------------------------------------------
 
void GmshReader::readMSH4(std::istream& is,
                           mesh::Mesh& mesh,
                           const std::unordered_map<int, Int>& pgMap)
{
    mesh::Data& d = mesh.data;
 
    // -----------------------------------------------------------------------
    // PhysicalNames (optional, used to resolve physical group → dim/tag)
    // -----------------------------------------------------------------------
 
    // Map: (dim, physTagGmsh) → user-facing name (for informational purposes)
    // We don't need the names for data output, but we do need the entity
    // physical-group tags to stamp RNG.
    // Nothing to do here beyond noting we'll get tags from $Entities.
 
    // -----------------------------------------------------------------------
    // Entities: build (dim, entityTag) → list-of-physicalGroupTags
    // -----------------------------------------------------------------------
 
    // entity physical groups: key = (dim << 32 | entityTag), value = physTags
    std::unordered_map<uint64_t, std::vector<int>> entityPhysTags;
 
    {
        // Seek to $Entities
        std::string line;
        {
            std::streampos saved = is.tellg();
            while (std::getline(is, line))
                if (line.find("$Entities") != std::string::npos) break;
            if (is.eof()) {
                // $Entities is optional (some Gmsh files omit it)
                is.clear();
                is.seekg(saved);
            } else {
                // Read entity counts: numPoints numCurves numSurfaces numVolumes
                Index np, nc, ns, nv;
                is >> np >> nc >> ns >> nv;
                std::getline(is, line); // consume rest
 
                auto readEntities = [&](Index count, int dim) {
                    for (Index i = 0; i < count; ++i) {
                        std::getline(is, line);
                        std::istringstream ss(line);
                        int tag;
                        ss >> tag;
                        // Skip bounding box: points have 3 coords, others have 6
                        if (dim == 0) {
                            double x, y, z; ss >> x >> y >> z;
                        } else {
                            double mn0, mn1, mn2, mx0, mx1, mx2;
                            ss >> mn0 >> mn1 >> mn2 >> mx0 >> mx1 >> mx2;
                        }
                        int numPhysTags; ss >> numPhysTags;
                        std::vector<int> ptags(numPhysTags);
                        for (int k = 0; k < numPhysTags; ++k) ss >> ptags[k];
                        uint64_t key = (uint64_t(dim) << 32) | uint64_t(tag);
                        entityPhysTags[key] = std::move(ptags);
                    }
                };
 
                readEntities(np, 0);
                readEntities(nc, 1);
                readEntities(ns, 2);
                readEntities(nv, 3);
            }
        }
    }
 
    // -----------------------------------------------------------------------
    // Nodes
    // -----------------------------------------------------------------------
 
    {
        is.clear(); is.seekg(0);
        std::string line;
        while (std::getline(is, line))
            if (line.find("$Nodes") != std::string::npos) break;
        if (is.eof())
            throw std::runtime_error("GmshReader MSH4: missing $Nodes section");
 
        // MSH4: numEntityBlocks numNodes minNodeTag maxNodeTag
        Index numEntityBlocks, numNodes, minTag, maxTag;
        is >> numEntityBlocks >> numNodes >> minTag >> maxTag;
 
        d.numNodes   = numNodes;
        d.spatialDim = 3;
        d.xyz.resize(numNodes * 3);
 
        // node tag → 0-based index
        std::unordered_map<Index, Index> nodeTagToIdx;
        nodeTagToIdx.reserve(numNodes);
 
        for (Index b = 0; b < numEntityBlocks; ++b) {
            int    entityDim, entityTag, parametric;
            Index  numNodesInBlock;
            is >> entityDim >> entityTag >> parametric >> numNodesInBlock;
 
            std::vector<Index> tags(numNodesInBlock);
            for (Index k = 0; k < numNodesInBlock; ++k)
                is >> tags[k];
 
            for (Index k = 0; k < numNodesInBlock; ++k) {
                double x, y, z;
                is >> x >> y >> z;
                Index idx = tags[k] - 1; // Use 0-based index derived from tag
                // For non-sequential tags build proper map below
                nodeTagToIdx[tags[k]] = 0; // placeholder; fix after all blocks
                d.xyz[k * 3 + 0] = static_cast<Real>(x); // temporary; re-indexed below
                d.xyz[k * 3 + 1] = static_cast<Real>(y);
                d.xyz[k * 3 + 2] = static_cast<Real>(z);
            }
        }
        // The above approach with non-sequential tags would misplace coords.
        // Re-read properly:
        is.clear(); is.seekg(0);
        while (std::getline(is, line))
            if (line.find("$Nodes") != std::string::npos) break;
 
        is >> numEntityBlocks >> numNodes >> minTag >> maxTag;
 
        std::vector<Real> xyzTmp(numNodes * 3);
        nodeTagToIdx.clear();
        Index nodeCounter = 0;
 
        for (Index b = 0; b < numEntityBlocks; ++b) {
            int   entityDim, entityTag, parametric;
            Index numNodesInBlock;
            is >> entityDim >> entityTag >> parametric >> numNodesInBlock;
 
            std::vector<Index> tags(numNodesInBlock);
            for (Index k = 0; k < numNodesInBlock; ++k)
                is >> tags[k];
 
            for (Index k = 0; k < numNodesInBlock; ++k) {
                double x, y, z;
                is >> x >> y >> z;
                Index idx = nodeCounter++;
                nodeTagToIdx[tags[k]] = idx;
                xyzTmp[idx * 3 + 0] = static_cast<Real>(x);
                xyzTmp[idx * 3 + 1] = static_cast<Real>(y);
                xyzTmp[idx * 3 + 2] = static_cast<Real>(z);
            }
        }
        d.xyz = std::move(xyzTmp);
 
        // Store nodeTagToIdx for use in Elements section
        // (passed via capture in the lambda below)
        // We'll re-use the variable in scope.
 
        // -----------------------------------------------------------------------
        // Elements
        // -----------------------------------------------------------------------
 
        {
            is.clear(); is.seekg(0);
            while (std::getline(is, line))
                if (line.find("$Elements") != std::string::npos) break;
            if (is.eof())
                throw std::runtime_error("GmshReader MSH4: missing $Elements section");
 
            Index numElemEntityBlocks, numElements, minElemTag, maxElemTag;
            is >> numElemEntityBlocks >> numElements >> minElemTag >> maxElemTag;
 
            struct RawElem {
                int          gmshType;
                int          physTag;   // resolved from entity
                std::vector<Index> nodeIdxs; // 0-based
            };
 
            std::vector<RawElem> volumeElems;
            std::vector<RawElem> boundaryElems;
            int detectedDim = 0;
 
            for (Index b = 0; b < numElemEntityBlocks; ++b) {
                int   entityDim, entityTag, gmshType;
                Index numElemsInBlock;
                is >> entityDim >> entityTag >> gmshType >> numElemsInBlock;
 
                const ElemTypeInfo* info = lookupElemType(gmshType);
 
                // Resolve physical tag for this entity
                int physTag = 0;
                uint64_t key = (uint64_t(entityDim) << 32) | uint64_t(entityTag);
                auto eit = entityPhysTags.find(key);
                if (eit != entityPhysTags.end() && !eit->second.empty())
                    physTag = eit->second[0]; // take first physical tag
 
                for (Index k = 0; k < numElemsInBlock; ++k) {
                    Index elemTag;
                    is >> elemTag;
 
                    if (!info) {
                        // Skip: consume node tags
                        // We don't know numNodes, so skip the line
                        std::getline(is, line);
                        continue;
                    }
 
                    std::vector<Index> nodeIdxs(info->numNodes);
                    for (Index ni = 0; ni < info->numNodes; ++ni) {
                        Index gmshNodeTag;
                        is >> gmshNodeTag;
                        auto nit = nodeTagToIdx.find(gmshNodeTag);
                        if (nit == nodeTagToIdx.end())
                            throw std::runtime_error(
                                "GmshReader MSH4: unknown node tag " +
                                std::to_string(gmshNodeTag));
                        nodeIdxs[ni] = nit->second;
                    }
 
                    if (info->dim >= 2) {
                        detectedDim = std::max(detectedDim, info->dim);
                        volumeElems.push_back({ gmshType, physTag, std::move(nodeIdxs) });
                    } else if (info->dim == 1) {
                        boundaryElems.push_back({ gmshType, physTag, std::move(nodeIdxs) });
                    }
                }
            }
 
            if (volumeElems.empty())
                throw std::runtime_error("GmshReader MSH4: no volume elements found");
 
            for (const auto& e : volumeElems)
                if (e.gmshType != volumeElems[0].gmshType)
                    throw std::runtime_error(
                        "GmshReader MSH4: mixed element types not supported");
 
            const ElemTypeInfo* typeInfo = lookupElemType(volumeElems[0].gmshType);
 
            d.parametricDim   = static_cast<Index>(detectedDim);
            d.spatialDim      = static_cast<Index>(detectedDim);
            d.numElements     = volumeElems.size();
            d.nodesPerElement = typeInfo->numNodes;
            d.facesPerElement = typeInfo->numFaces;
            d.basisOrder      = std::vector<Index>(detectedDim, 1);
 
            // Trim xyz to spatialDim
            if (d.spatialDim == 2) {
                std::vector<Real> xyz2(d.numNodes * 2);
                for (Index n = 0; n < d.numNodes; ++n) {
                    xyz2[n * 2 + 0] = d.xyz[n * 3 + 0];
                    xyz2[n * 2 + 1] = d.xyz[n * 3 + 1];
                }
                d.xyz = std::move(xyz2);
            }
 
            // IEN
            d.ien.resize(d.numElements * d.nodesPerElement);
            for (Index e = 0; e < d.numElements; ++e)
                for (Index k = 0; k < d.nodesPerElement; ++k)
                    d.ien[e * d.nodesPerElement + k] = volumeElems[e].nodeIdxs[k];
 
            // RNG
            d.rng.assign(d.numElements * d.facesPerElement, -1);
 
            // Build face map (same logic as MSH2)
            std::vector<std::vector<Index>> faceDefs;
            if (d.nodesPerElement == 4)
                faceDefs = { {0,1}, {1,2}, {2,3}, {3,0} };
            else if (d.nodesPerElement == 8)
                faceDefs = { {0,1,2,3}, {4,5,6,7}, {0,1,5,4},
                             {1,2,6,5}, {2,3,7,6}, {3,0,4,7} };
 
            using NodeSet = std::vector<Index>;
            auto nodeSetToKey = [](NodeSet s) -> std::string {
                std::sort(s.begin(), s.end());
                std::string k;
                k.reserve(s.size() * sizeof(Index));
                for (Index v : s)
                    k.append(reinterpret_cast<const char*>(&v), sizeof(Index));
                return k;
            };
 
            struct FaceRef { Index elem; Index face; };
            std::unordered_map<std::string, FaceRef> faceMap;
            faceMap.reserve(d.numElements * d.facesPerElement);
 
            for (Index e = 0; e < d.numElements; ++e) {
                for (Index f = 0; f < d.facesPerElement; ++f) {
                    NodeSet ns;
                    for (Index li : faceDefs[f])
                        ns.push_back(d.ien[e * d.nodesPerElement + li]);
                    faceMap[nodeSetToKey(ns)] = { e, f };
                }
            }
 
            for (const auto& be : boundaryElems) {
                NodeSet ns(be.nodeIdxs.begin(), be.nodeIdxs.end());
                auto key = nodeSetToKey(ns);
                auto fit = faceMap.find(key);
                if (fit == faceMap.end()) continue;
 
                Int solverTag = be.physTag;
                if (!pgMap.empty()) {
                    auto pit = pgMap.find(be.physTag);
                    if (pit != pgMap.end()) solverTag = pit->second;
                }
                d.rng[fit->second.elem * d.facesPerElement + fit->second.face] = solverTag;
            }
        }
    }
}
