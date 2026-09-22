# Architecture

<!-- STUB — outline below, full prose to be written. -->

A narrative walkthrough of the codebase, aimed at a human reading top-to-bottom rather than
looking something up. See `CLAUDE.md` for the terse structural reference this expands on.

## Planned sections

1. **The layering, and why it's bottom-up.** `mesh`/`io` → `linalg`/`fem` → `equation` →
   `solver` → `application`. What each layer is allowed to know about the ones above it (nothing).
2. **A traced example run.** Follow `examples/heateq/steady/constant_conductivity/2d/config.yaml` from
   `HeatApplication::main()` through `HeatConfigParser` → `HeatDispatcher` (runtime config
   resolved into a concrete template instantiation) → `solver::stage::SteadyStage<HeatProblemT>` →
   `fem::assembly::Assembler` → `linalg::solver::iterative::cg::Solver` → VTK output. This is the
   single best orientation for a new reader.
3. **Why concepts, not virtual dispatch.** The cost/benefit of the C++20 `concept`-constrained
   extension-point design (`fem/form`, `fem/evaluator`, `solver/stage`, `linalg/operator`) versus
   an OOP interface hierarchy — compile-time monomorphization for the hot assembly loop, at the
   cost of longer build times and template-error surface area (mitigated by the concept
   constraints themselves turning most of those into a one-line diagnostic).
4. **Backend-agnostic design.** How the `Backend` tag pattern (`linalg/types/backend`,
   `fem/assembly/backend`, `fem/boundary/backend`) is meant to make a CUDA specialization a
   drop-in — and what "drop-in" actually requires of new code (see `CLAUDE.md` → Extensibility
   priorities).
5. **Runtime dispatch vs. compile-time templates.** How `fem::dispatch` collapses the
   (basis, quadrature, backend) config space into a bounded set of template instantiations
   resolved once at startup, so the assembly loop itself stays monomorphic.
6. **`SegregatedStage` and multi-physics composition.** How a composite stage lets `Steady`/
   `Transient` stay ignorant of whether they're driving one physics or a segregated sequence of
   several (e.g. a future pressure-velocity + turbulence RANS coupling).
