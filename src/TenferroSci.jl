"""
    TenferroSci

Special functions from [tenferro-sci-rs](https://github.com/tensor4all/tenferro-sci-rs),
called through [RustCall.jl](https://github.com/AtelierArith/RustCall.jl).

The Rust half is `rust/tenferro_sci_ffi`, a thin `#[julia]` wrapper around the
upstream `specials` crate. `specials` is generic over a `Real` scalar with two
evaluation strategies:

* `Float64` uses kernels ported from SciPy's xsf, so the results agree with
  `scipy.special` to within a few ulp.
* The double-double path (about 32 significant digits) uses
  precision-adaptive algorithms written from DLMF. Pass `Float64x2` from
  MultiFloats.jl (re-exported here) — converted to `Df64` internally, with the
  result converted back — or `TenferroSci.Df64` directly.

```julia
julia> gamma(5.0)
24.0

julia> zeta(2.0, 1.0) ≈ pi^2 / 6
true

julia> y = gamma(Float64x2(2.5)); big(y)
1.3293403881791370204736256125058588870981620920917903461603558...

julia> isapprox(big(y), big"1.3293403881791370204736256125058588870981620920917903461603558"; rtol=1e-30)
true
```

Domain errors and overflow are values, not exceptions: the result is `NaN`
outside the domain and at poles whose sign is undefined, and `±Inf` at poles
with a definite sign or on overflow, as in SciPy.

The `Float64` methods accept `Float64`, `Float32`, `Float16`, integers,
rationals and irrationals. `BigFloat` is refused rather than silently rounded
to 53 bits; use `Df64` for about 32 digits, or `SpecialFunctions.jl` when you
need `BigFloat`.

The exported names are `scipy.special` names. `Base` exports none of them
except through `SpecialFunctions.jl`, so qualify the call when both packages
are loaded.

# Scope

This release binds the gamma family of `specials`:

`gamma`, `lgamma`, `gammasgn`, `rgamma`, `digamma`, `zeta`, `beta`, `lbeta`,
`poch`, `lgam1p`, `binom`, `zetac`, `trigamma`, `digammainv`, `gen_harmonic`,
`stirling2`.

The other families re-exported by `specials` (`fermi_dirac`, `bose_einstein`,
`debye`, `coulomb`, …) are not bound yet.
"""
module TenferroSci

using MultiFloats: Float64x2
using RustCall

# The crate is scanned, built and bound while this package is precompiled.
# `submodule` names the generated module `TenferroSci.Bindings`; the bindings
# are implementation detail, the public API is defined below.
@rust_crate joinpath(@__DIR__, "..", "rust", "tenferro_sci_ffi") submodule = "Bindings"

export gamma, lgamma, gammasgn, rgamma, digamma, zeta, beta, lbeta, poch,
    lgam1p, binom, zetac, trigamma, digammainv, gen_harmonic, stirling2,
    Float64x2

include("df64.jl")
include("functions.jl")
include("multifloats.jl")

end # module TenferroSci
