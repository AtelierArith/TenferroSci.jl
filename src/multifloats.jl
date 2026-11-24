# `Float64x2` is the user-facing high-precision type. Every method here
# converts its `Float64x2` arguments to `Df64`, calls the same double-double
# method the rest of the package uses, and converts the result back.
# `Float64x2` is exactly the two-limb representation `hi + lo`, so both
# conversions are exact and the result keeps the ~32 digits of the Rust
# algorithms.
#
# The Rust algorithms are two-limb, so only `Float64x2` is accepted. Rounding a
# `Float64x3` or `Float64x4` to `Float64x2` is explicit: call
# `Float64x2(x)` before the special function, and the result has about 32
# digits, not 48 or 64.

const MF2 = Float64x2

"""
    Df64(x::Float64x2)

Exact conversion. `Float64x2` values are normalised, so this only copies the
two limbs.
"""
Df64(x::MF2) = Df64(getfield(x, :_limbs)...)

"""
    Float64x2(x::Df64)

Exact conversion. The two limbs of `x` become the two limbs of the result; a
non-finite `Df64` keeps its `Inf`/`NaN` (MultiFloats itself turns `Inf` into
`NaN` only in arithmetic).
"""
MF2(x::Df64) = MF2((x.hi, x.lo))

const UNARY = (
    :gamma, :lgamma, :gammasgn, :rgamma, :digamma, :lgam1p, :zetac, :trigamma,
    :digammainv,
)

const BINARY = (
    :zeta, :beta, :lbeta, :poch, :binom, :gen_harmonic, :stirling2,
)

for f in UNARY
    @eval TenferroSci.$f(x::MF2) = MF2(TenferroSci.$f(Df64(x)))
end

for f in BINARY
    @eval begin
        TenferroSci.$f(a::MF2, b::MF2) = MF2(TenferroSci.$f(Df64(a), Df64(b)))
        TenferroSci.$f(a::MF2, b::F64Input) = MF2(TenferroSci.$f(Df64(a), Float64(b)))
        TenferroSci.$f(a::F64Input, b::MF2) = MF2(TenferroSci.$f(Float64(a), Df64(b)))
        TenferroSci.$f(a::MF2, b::Df64) = MF2(TenferroSci.$f(Df64(a), b))
        TenferroSci.$f(a::Df64, b::MF2) = MF2(TenferroSci.$f(a, Df64(b)))
    end
end
