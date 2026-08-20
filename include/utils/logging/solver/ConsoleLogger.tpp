namespace pdesolver::utils::logging {

	// log()
	template<typename DataType>
	void ConsoleLogger::log(Index iter, const std::vector<DataType>& perDOFAbs, DataType flopsThisIter) const {

		// capture the relative-residual baseline and start the clock on the first call
		if (!baselineCaptured_) {
			initialPerDOFNorms_.assign(perDOFAbs.begin(), perDOFAbs.end());
			baselineCaptured_ = true;
			startTime_ = std::chrono::steady_clock::now();
		}

		// lifetime accounting, used only by summary() -- independent of print interval
		totalFlops_ += static_cast<double>(flopsThisIter);
		lastIter_ = iter;

		// convert to relative, per-DOF, against the iteration-0 baseline
		lastRelPerDOF_.assign(perDOFAbs.size(), Real(0));
		for (Index i = 0; i < perDOFAbs.size(); ++i) {
			Real base = initialPerDOFNorms_[i];
			lastRelPerDOF_[i] = (base > Real(0)) ? static_cast<Real>(perDOFAbs[i]) / base : Real(0);
		}

		// check to print based on config
		if (interval == 0) return;
		if ((iter > 0) && ((iter % interval) != 0)) return;

		// print header
		if (printHeader && iter == 0) {
			printBanner();
			printColumnHeader();
			const_cast<ConsoleLogger*>(this)->printHeader = false;
		}

		const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime_).count();

		// set print format
		std::cout << "  " << std::left << std::setw(6) << iter << "  ";
		std::cout << std::scientific << std::setprecision(4);

		for (const auto& v : lastRelPerDOF_) {
			std::cout << std::setw(14) << v << "  ";
		}

		std::cout << std::fixed << std::setprecision(3) << elapsed << "s";
		std::cout << "\n";

	}

	// computePerDOFNorms() helper
	template<typename DataType>
	std::vector<DataType> ConsoleLogger::computePerDOFNorms(const DataType* r, Index totalSize) const {

		const Index D = dofNames.size();
		std::vector<DataType> norms(D > 0 ? D : 1, DataType(0));

		// single field: the whole residual belongs to the one DOF
		if (D <= 1 || freeDOFsPerField == 0) {
			DataType sum = DataType(0);
			for (Index i = 0; i < totalSize; ++i) sum += r[i] * r[i];
			norms[0] = std::sqrt(sum);
			return norms;
		}

		// split total residual norm to each dof component
		if (dofOrdering == fem::dof::DOFOrdering::Interleaved) {

			// r[i] belongs to component (i % D)
			for (Index i = 0; i < totalSize; ++i) {
				Index comp = (i % D);
				norms[comp] += r[i] * r[i];
			}

		} else {

			// r[c*freeDOFsPerField .. (c+1)*freeDOFsPerField]
			for (Index c = 0; c < D; ++c){
				Index start = c * freeDOFsPerField;
				Index end = start + freeDOFsPerField;
				if (end > totalSize) end = totalSize;
				for (Index i = start; i < end; ++i){
					norms[c] += r[i] * r[i];
				}
			}

		}

		// take 2-norm of each component
		for (auto& v : norms) v = std::sqrt(v);

		return norms;

	}

	inline void ConsoleLogger::summary(bool converged) const {

		const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime_).count();
		const double avgFlopsPerIter = (lastIter_ > 0) ? (totalFlops_ / static_cast<double>(lastIter_)) : 0.0;

		std::cout << "  [" << solverName << "] " << (converged ? "converged" : "did not converge") << ": "
				  << lastIter_ << " iterations, ";

		for (Index i = 0; i < lastRelPerDOF_.size(); ++i) {
			std::cout << "res[" << dofNames[i] << "] " << std::scientific << std::setprecision(4) << lastRelPerDOF_[i];
			if (i + 1 < lastRelPerDOF_.size()) std::cout << ", ";
		}

		std::cout << ", avg " << std::scientific << std::setprecision(2) << avgFlopsPerIter << " flops/iter, "
				  << std::fixed << std::setprecision(3) << elapsed << "s\n";

	}

	inline void ConsoleLogger::printBanner() const {

		const int width = 60;

		std::cout << "\n";
		std::cout << "  " << std::string(width, '=') << "\n";
		std::cout << "  " << equationName << " - " << solverName << "\n";
		std::cout << "  Preconditioner: " << preconditionerName << "\n";

		for (const auto& param : extraParams) {
			std::cout << "  " << param.first << ": " << param.second << "\n";
		}

		if (dofNames.size() > 1) {
			std::cout << "  Components:";
			for (const auto& n : dofNames) std::cout << "  " << n;
			std::cout << "\n";
			std::cout << "  DOF ordering: "
					  << (dofOrdering == fem::dof::DOFOrdering::Interleaved
						  ? "Interleaved (node-major)"
						  : "Block (field-major)")
					  << "\n";
		}

		std::cout << "  " << std::string(width, '=') << "\n";

	}

	inline void ConsoleLogger::printColumnHeader() const {

		std::cout << "  " << std::left << std::setw(6) << "Iter" << "  ";

		for (const auto& name : dofNames) {
			std::cout << std::setw(14) << ("Res[" + name + "]") << "  ";
		}

		std::cout << std::setw(10) << "Elapsed" << "\n";

		const Index width = 8 + 16 * dofNames.size() + 10;
		std::cout << "  " << std::string(std::max<Index>(30, width), '-') << "\n";

	}

} // namespace pdesolver::utils::logging
