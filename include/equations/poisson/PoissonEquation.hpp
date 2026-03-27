#ifndef POISSON_EQUATION
#define POISSON_EQUATION

#include "equations/poisson/boundary/BoundaryFluxFunction.hpp"
#include "equations/poisson/boundary/BoundaryValueFunction.hpp"

#include "equations/poisson/eval/SourceFunction.hpp"
#include "equations/poisson/eval/DefaultModel.hpp"
#include "equations/poisson/eval/ConductivityModel.hpp"
#include "equations/poisson/eval/EvalElement.hpp"
#include "equations/poisson/eval/EvalField.hpp"
#include "equations/poisson/eval/EvalQuadraturePointVolume.hpp"
#include "equations/poisson/eval/EvalQuadraturePointBoundary.hpp"

#include "equations/poisson/form/DiffusionForm.hpp"
#include "equations/poisson/form/SourceForm.hpp"
#include "equations/poisson/form/FluxBoundaryForm.hpp"

#endif
