# TenferroSci.jl

Julia bindings for the special functions of
[tenferro-sci-rs](https://github.com/tensor4all/tenferro-sci-rs), called through
[RustCall.jl](https://github.com/AtelierArith/RustCall.jl).

The Rust side is [`specials`](https://github.com/tensor4all/tenferro-sci-rs/tree/main/crates/specials),
the tenferro-free core of tenferro-sci-rs. It evaluates every function in two
ways:

* `Float64` uses kernels ported from SciPy's xsf, so results agree with
  `scipy.special` to within a few ulp.
* `Float64x2` (from MultiFloats.jl, about 32 significant digits) selects the
  double-double path, which uses precision-adaptive algorithms written from
  DLMF. Internally the value is the two-limb `TenferroSci.Df64`; the upstream
  test suite checks that path against mpmath with tolerances around 1e-29.

This package exposes both paths with the same names as `scipy.special`.

## Status

The gamma family of `specials` is bound:

`gamma`, `lgamma`, `gammasgn`, `rgamma`, `digamma`, `zeta`, `beta`, `lbeta`,
`poch`, `lgam1p`, `binom`, `zetac`, `trigamma`, `digammainv`, `gen_harmonic`,
`stirling2`.

The other families re-exported by `specials` (`fermi_dirac`, `bose_einstein`,
`debye`, `coulomb`, `whittaker`, …) are not bound yet. See
[Roadmap](#Roadmap).

The `specials` crate is not on crates.io, so the Rust crate pins the upstream
commit in `rust/tenferro_sci_ffi/Cargo.toml`.

## Installation

Julia 1.12 or later is required (RustCall.jl uses Julia's `@rust_crate`
precompilation support). A Rust toolchain is not required on `PATH`:
RustCall.jl installs one through Julia artifacts when needed.

```julia
using Pkg
Pkg.develop(url = "https://github.com/atelierarith/TenferroSci.jl")
```

The first load scans and compiles the Rust crate. Later loads use RustCall's
cache; editing the Rust sources triggers a rebuild.

## Usage

```julia
julia> using TenferroSci

julia> gamma(5.0)
24.0

julia> zeta(2.0, 1.0) ≈ pi^2 / 6
true

julia> beta(2.0, 3.0)
0.08333333333333333

julia> stirling2(10.0, 7.0)
5880.0
```

### High precision

The high-precision interface is `Float64x2` from
[MultiFloats.jl](https://github.com/dzhang314/MultiFloats.jl), re-exported by
this package. The argument is converted to the internal `Df64` type, the
double-double method of the special function runs, and the result is converted
back to `Float64x2`. Both conversions are exact — the two types are the same
two-limb pair `hi + lo` — so the result keeps the ~32 digits of the Rust
algorithms and supports MultiFloats' correctly rounded arithmetic:

```julia
julia> using TenferroSci

julia> x = Float64x2(2.5);

julia> y = gamma(x);

julia> y isa Float64x2
true

julia> y * y / x
3.5000000000000000000000000000000000

julia> isapprox(BigFloat(y), big"1.3293403881791370204736256125058588870981620920917903461603558"; rtol = 1e-30)
true
```

`TenferroSci.Df64`, `Float64x2`, `Float64` and integers can be mixed in binary
calls; the result is a `Float64x2` when either argument is:

```julia
julia> zeta(Float64x2(2), 1)
Float64x2(1.6449340668482264)
```

Only `Float64x2` is accepted. The Rust algorithms are two-limb, so a
`Float64x3` or `Float64x4` must be rounded explicitly with `Float64x2(x)`, and
the result then has about 32 digits, not 48 or 64. Conversion preserves `Inf`
and `NaN`; MultiFloats itself turns `Inf` into `NaN` in arithmetic, so do not
continue arithmetic from a special value.

#### Internal representation: `TenferroSci.Df64`

`Df64` is the two-limb type used at the Rust boundary. It is not exported, but
remains available as `TenferroSci.Df64` for callers that want the raw
representation:

```julia
julia> x = TenferroSci.Df64(2.5);

julia> y = gamma(x);

julia> big(y)
1.329340388179137020473625612505858887098162092091790346160355...

julia> isapprox(big(y), big"1.3293403881791370204736256125058588870981620920917903461603558"; rtol = 1e-30)
true
```

`Df64` is not a subtype of `Real` and has no arithmetic of its own. It supports
the functions here, conversion with `big`, `Float64` and `Float64x2`, and
`isnan`, `isinf`, `isfinite`, `iszero`, `signbit`. It is an immutable pair
`hi + lo`; the constructor normalises any finite pair with a two-sum.

The `Float64` methods accept `Float64`, `Float32`, `Float16`, integers,
rationals and irrationals. `BigFloat` is refused rather than silently rounded
to 53 bits; use `Float64x2`/`Df64` for about 32 digits, or
`SpecialFunctions.jl` when you need `BigFloat`.

Domain errors and overflow are values, not exceptions, following SciPy: `NaN`
outside the domain and at poles whose sign is undefined, `±Inf` at poles with a
definite sign and on overflow.

When `SpecialFunctions.jl` is also loaded, qualify calls: both packages export
the same `scipy.special` names.

## How it works

```
Julia: TenferroSci.gamma(x)
        │  Float64 → f64 method
        │  Float64x2 → Df64 → *_dd method → Float64x2
        ▼
Rust:  rust/tenferro_sci_ffi/src/lib.rs   (#[julia] wrappers, concrete types)
        │
        ▼
       specials::gamma::<f64>  or  specials::gamma::<xprec::Df64>
```

`specials` is generic over `Real`, but `#[julia]` cannot bind a generic Rust
function: an `extern "C"` entry point needs concrete types. The wrapper crate
therefore instantiates each function twice, once for `f64` and once for
`xprec::Df64`.

RustCall's FFI contract does not return a `#[julia]` struct by value, so the
double-double functions take a `&mut Dd` out-parameter. The Julia wrapper
allocates one `Dd` handle per call, reads `hi` and `lo`, and lets the handle be
collected. Each `Df64` call therefore performs one Rust allocation; the
arithmetic itself dominates.

The generated low-level module is available as `TenferroSci.Bindings` — for
example `TenferroSci.Bindings.gamma(2.5)` calls the raw `f64` binding without
the Julia conversion layer. It is implementation detail and may change.

## Testing

```julia
using Pkg
Pkg.test()
```

`test/reference_values.csv` holds mpmath values at 120 decimal digits for 79
argument sets (`uv run test/gen_reference.py` regenerates it). The test suite
checks the `Float64` path against them with xsf-level tolerances and the `Df64`
path with 5e-27 relative tolerance, plus the special values, argument
conversion, and the `Df64` constructor.

These tests cover the FFI boundary (argument order, `hi`/`lo` packing, return
values). The algorithms themselves are validated upstream in tenferro-sci-rs
against SciPy's xsref tables and mpmath; see the `specials` crate documentation.

Rust side:

```bash
cd rust/tenferro_sci_ffi
cargo fmt --check
cargo clippy --release -- -D warnings
cargo build --release
```

## Roadmap

* Bind the remaining `specials` families re-exported from
  `few_special_functions` (Fermi–Dirac, Bose–Einstein, Debye, Coulomb,
  Whittaker, Marcum Q, …).
* Follow upstream as it adds error functions, Bessel functions, incomplete
  gamma/beta, elliptic integrals and the rest of the SciPy special-function
  set.
* Support `Float64x3`/`Float64x4` natively. The Rust algorithms truncate at the
  scalar type's epsilon, so this needs a three- or four-limb Rust type
  implementing `specials::Real`; until then, round to `Float64x2` explicitly.
* Once tenferro-sci-rs publishes `tenferro-krylov`, `tenferro-optimize` and the
  other crates, bind them the same way.

## License

The Julia package and the FFI crate are MIT OR Apache-2.0. The ported code in
`specials` keeps its upstream licenses (xsf/Cephes BSD-3-Clause, BSL-1.0, MIT);
see the upstream crate.
