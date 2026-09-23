# Architecture

<!-- Outline — topics to discuss/draft, not final content. -->

- Layering and why it's bottom-up (`mesh`/`io` → `linalg`/`fem` → `equation` → `solver` →
  `application`)
- A traced example run: config → `Dispatcher` → templated `Stage` → `Assembler` → linear solver →
  output
- Why concepts over virtual dispatch (monomorphic hot loop vs. build time / error surface)
- Backend-agnostic design: the `Backend` tag pattern, what a CUDA drop-in actually requires
- Runtime dispatch vs. compile-time templates (`fem::dispatch`)
- `SegregatedStage` and multi-physics composition
