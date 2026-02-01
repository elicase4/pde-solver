#ifndef PDESOLVER_FEM
#define PDESOLVER_FEM

#include "fem/assmebly/Assembler.hpp"

#include "fem/boundary/BoundaryApplicator.hpp"
#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/boundary/BoundaryRegistry.hpp"

#include "fem/dof/AlgebraicDOF.hpp"

#include "fem/eval/EvalElement.hpp"
#include "fem/eval/EvalField.hpp"

#include "fem/form/BilinearForm.hpp"
#include "fem/form/LinearForm.hpp"
#include "fem/form/NonlinearForm.hpp"
#include "fem/form/NonlinearTangentForm.hpp"

#include "fem/geometry/JacobianTransform.hpp"

#include "fem/basis/Lagrange1D.hpp"
#include "fem/basis/LagrangeQuad.hpp"
#include "fem/basis/LagrangeHex.hpp"

#include "fem/quadrature/GaussQuadrature1D.hpp"
#include "fem/quadrature/GaussQuadratureQuad.hpp"
#include "fem/quadrature/GaussQuadratureHex.hpp"

#endif
