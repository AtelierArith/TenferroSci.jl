# Public API. Each function has an `f64` method over `Real` (the xsf port)
# and a `Df64` method over the double-double path. The generated
# `TenferroSci.Bindings.<name>_dd` writes its result into a Rust-owned `Dd`
# handle, which is unpacked immediately; the handle is collected by the GC.

# Accepted by the `Float64` methods. `BigFloat` is deliberately absent: it
# would be silently rounded to 53 bits. Use `Df64` for about 32 digits.
const F64Input = Union{Float16,Float32,Float64,Integer,Rational,Irrational}

@inline _dd1(f, x::Df64) = unpack_dd(f(x.hi, x.lo))

@inline _dd2(f, x::Df64, y::Df64) = unpack_dd(f(x.hi, x.lo, y.hi, y.lo))

"""
    gamma(x)

Γ(x), the gamma function (DLMF 5.2.1), like `scipy.special.gamma`.

Poles: `gamma(±0) = ±Inf`; at negative integers and at `-Inf` the result is
`NaN`. `gamma(Inf) = Inf`, and the result overflows to `Inf` above
x ≈ 171.62.

The `Float64` method calls the xsf port and the `Df64` method the
precision-adaptive one.

# Examples

```julia
julia> gamma(5.0)
24.0

julia> isnan(gamma(-1.0))
true

julia> isapprox(big(gamma(Df64(2.5))), gamma(big"2.5"); rtol=1e-30)
true
```
"""
gamma(x::F64Input) = Bindings.gamma(Float64(x))
gamma(x::Df64) = _dd1(Bindings.gamma_dd, x)

"""
    lgamma(x)

ln|Γ(x)|, like `scipy.special.gammaln`.

`Inf` at the poles (zero and the negative integers). Near its zeros at x = 1
and x = 2 the `Float64` kernel is accurate in absolute rather than relative
terms (about 1e-16); the `Df64` method stays relatively accurate there.
"""
lgamma(x::F64Input) = Bindings.lgamma(Float64(x))
lgamma(x::Df64) = _dd1(Bindings.lgamma_dd, x)

"""
    gammasgn(x)

The sign of Γ(x): `1.0` or `-1.0`, `±1` at `±0`, `NaN` at the poles and at
`-Inf`, like `scipy.special.gammasgn`.
"""
gammasgn(x::F64Input) = Bindings.gammasgn(Float64(x))
gammasgn(x::Df64) = _dd1(Bindings.gammasgn_dd, x)

"""
    rgamma(x)

1/Γ(x), the reciprocal gamma function, like `scipy.special.rgamma`. Entire,
with zeros at the poles of `gamma`.

Results below `floatmin(Float64)` (x ≳ 171.6) flush to zero.
"""
rgamma(x::F64Input) = Bindings.rgamma(Float64(x))
rgamma(x::Df64) = _dd1(Bindings.rgamma_dd, x)

"""
    digamma(x)

ψ(x) = Γ'(x)/Γ(x), the digamma function (DLMF 5.2.2), like
`scipy.special.digamma` and `scipy.special.psi`.

`digamma(±0) = ∓Inf`; `NaN` at the negative integers and at `-Inf`.
"""
digamma(x::F64Input) = Bindings.digamma(Float64(x))
digamma(x::Df64) = _dd1(Bindings.digamma_dd, x)

"""
    zeta(s, q)

ζ(s, q) = Σ_{k≥0} (k + q)^-s, the Hurwitz zeta function (DLMF 25.11.1), like
`scipy.special.zeta`.

For `q = 1` this is the Riemann zeta function for every real `s`, with the
analytic continuation for `s < 1`. For `q != 1` it is defined for `s > 1`:
`s = 1` gives `Inf`, `s < 1` gives `NaN`, a non-positive integer `q` gives
`Inf`, and `q < 0` requires an integer `s` (`NaN` otherwise).
"""
zeta(s::F64Input, q::F64Input) = Bindings.zeta(Float64(s), Float64(q))
zeta(s::Df64, q::F64Input) = zeta(s, Df64(q))
zeta(s::F64Input, q::Df64) = zeta(Df64(s), q)
zeta(s::Df64, q::Df64) = _dd2(Bindings.zeta_dd, s, q)

"""
    beta(a, b)

B(a, b) = Γ(a)Γ(b)/Γ(a+b), the beta function (DLMF 5.12.1), like
`scipy.special.beta`.

For a negative integer argument the result follows SciPy: finite when the
other argument is an integer that makes B a polynomial value, `Inf` otherwise.
"""
beta(a::F64Input, b::F64Input) = Bindings.beta(Float64(a), Float64(b))
beta(a::Df64, b::F64Input) = beta(a, Df64(b))
beta(a::F64Input, b::Df64) = beta(Df64(a), b)
beta(a::Df64, b::Df64) = _dd2(Bindings.beta_dd, a, b)

"""
    lbeta(a, b)

ln|B(a, b)|, like `scipy.special.betaln`. Accurate in absolute terms near the
zeros of ln|B|, where |B| = 1.
"""
lbeta(a::F64Input, b::F64Input) = Bindings.lbeta(Float64(a), Float64(b))
lbeta(a::Df64, b::F64Input) = lbeta(a, Df64(b))
lbeta(a::F64Input, b::Df64) = lbeta(Df64(a), b)
lbeta(a::Df64, b::Df64) = _dd2(Bindings.lbeta_dd, a, b)

"""
    poch(a, m)

(a)ₘ = Γ(a + m)/Γ(a), the Pochhammer symbol (DLMF 5.2.5), like
`scipy.special.poch`.

For a non-negative integer `m` this is the product `a (a+1) ⋯ (a+m-1)`. It is
`0` when `a` is a non-positive integer and `a + m` is not, and `Inf` the other
way round.
"""
poch(a::F64Input, m::F64Input) = Bindings.poch(Float64(a), Float64(m))
poch(a::Df64, m::F64Input) = poch(a, Df64(m))
poch(a::F64Input, m::Df64) = poch(Df64(a), m)
poch(a::Df64, m::Df64) = _dd2(Bindings.poch_dd, a, m)

"""
    lgam1p(x)

ln Γ(1 + x), accurate for `x` near 0 and near 1 where `lgamma(1 + x)` is not,
like `scipy.special.lgam1p`.
"""
lgam1p(x::F64Input) = Bindings.lgam1p(Float64(x))
lgam1p(x::Df64) = _dd1(Bindings.lgam1p_dd, x)

"""
    binom(n, k)

C(n, k), the binomial coefficient, extended to real arguments as
Γ(n+1)/(Γ(k+1)Γ(n-k+1)), like `scipy.special.binom`.

The `Float64` method carries xsf's error of up to about 3e-10 when
`k > 1e8 |n|` and `k < 1e6` (for example `binom(1e-9, 0.5)`); the `Df64`
method is accurate there.
"""
binom(n::F64Input, k::F64Input) = Bindings.binom(Float64(n), Float64(k))
binom(n::Df64, k::F64Input) = binom(n, Df64(k))
binom(n::F64Input, k::Df64) = binom(Df64(n), k)
binom(n::Df64, k::Df64) = _dd2(Bindings.binom_dd, n, k)

"""
    zetac(x)

ζ(x) − 1, like `scipy.special.zetac`. Accurate for large `x` where `zeta(x)`
rounds to 1.

`zetac(1) = Inf`, `zetac(Inf) = 0`, `NaN` at `-Inf`.
"""
zetac(x::F64Input) = Bindings.zetac(Float64(x))
zetac(x::Df64) = _dd1(Bindings.zetac_dd, x)

"""
    trigamma(x)

ψ'(x), the trigamma function (DLMF 5.15.1), like
`scipy.special.polygamma(1, x)`.

`Inf` at zero and the negative integers, `NaN` at `-Inf`. For `x > 0` it equals
`zeta(2, x)`.
"""
trigamma(x::F64Input) = Bindings.trigamma(Float64(x))
trigamma(x::Df64) = _dd1(Bindings.trigamma_dd, x)

"""
    digammainv(y)

The `x > 0` with `digamma(x) = y`, the inverse of [`digamma`](@ref) on its
increasing branch, like `scipy.special.digammainv`.

`0` for `y = -Inf`; `Inf` for `y` above ln(`floatmax(Float64)`) ≈ 709.78 and
for `y = Inf`.
"""
digammainv(y::F64Input) = Bindings.digammainv(Float64(y))
digammainv(y::Df64) = _dd1(Bindings.digammainv_dd, y)

"""
    gen_harmonic(n, a)

H(n, a) = Σ_{i=1}^{n} i^-a, the generalized harmonic number, like xsf's
`gen_harmonic`.

`n` must be a positive integer; a non-integer or `n < 1` gives `NaN`. `a` may
be any real number. `n = Inf` gives the limit: `zeta(a)` for `a > 1` and `Inf`
for `a <= 1`.
"""
gen_harmonic(n::F64Input, a::F64Input) = Bindings.gen_harmonic(Float64(n), Float64(a))
gen_harmonic(n::Df64, a::F64Input) = gen_harmonic(n, Df64(a))
gen_harmonic(n::F64Input, a::Df64) = gen_harmonic(Df64(n), a)
gen_harmonic(n::Df64, a::Df64) = _dd2(Bindings.gen_harmonic_dd, n, a)

"""
    stirling2(n, k)

S(n, k), the Stirling number of the second kind: the number of ways to
partition `n` elements into `k` non-empty subsets, like
`scipy.special.stirling2`.

`n` and `k` are truncated towards zero; `S(n, n) = 1`, `S(n, k) = 0` for
`k > n` and for `k = 0 < n`. For `n > 219` the `Float64` method is xsf's Temme
asymptotic approximation (relative error up to about 1e-6 at `n = 220`); the
`Df64` method evaluates the exact recurrence within its work limit.
"""
stirling2(n::F64Input, k::F64Input) = Bindings.stirling2(Float64(n), Float64(k))
stirling2(n::Df64, k::F64Input) = stirling2(n, Df64(k))
stirling2(n::F64Input, k::Df64) = stirling2(Df64(n), k)
stirling2(n::Df64, k::Df64) = _dd2(Bindings.stirling2_dd, n, k)
