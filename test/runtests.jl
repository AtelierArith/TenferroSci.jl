using Test
using TenferroSci
using TenferroSci: Df64

# Reference values: `test/gen_reference.py`, mpmath at 120 decimal digits.
# The arguments are the shortest round-tripping `Float64` representations, so
# Julia parses them back to exactly the values the Rust side receives.
struct ReferenceRow
    name::String
    a::Float64
    b::Union{Nothing,Float64}
    value::BigFloat
end

function read_reference()
    rows = ReferenceRow[]
    for line in eachline(joinpath(@__DIR__, "reference_values.csv"))
        startswith(line, '#') && continue
        isempty(strip(line)) && continue
        f = split(line, ',')
        name = String(f[1])
        a = parse(Float64, f[2])
        b = isempty(f[3]) ? nothing : parse(Float64, f[3])
        value = parse(BigFloat, f[4])
        push!(rows, ReferenceRow(name, a, b, value))
    end
    return rows
end

at_big_precision(f) = setprecision(f, BigFloat, 512)

const ROWS = at_big_precision(read_reference)

# `Float64` parity with `scipy.special`, through the xsf port. The tolerances
# are loose enough for the few documented xsf behaviours (binom's large-k
# branch, stirling2's Temme branch) and tight enough to catch an argument or
# result mix-up at the FFI boundary.
const F64_TOL = Dict(
    "gamma" => 5e-14,
    "lgamma" => 5e-14,
    "gammasgn" => 0.0,          # exact
    "rgamma" => 5e-14,
    "digamma" => 5e-14,
    "zeta" => 2e-13,
    "beta" => 2e-13,
    "lbeta" => 2e-13,
    "poch" => 1e-11,
    "lgam1p" => 1e-13,
    "binom" => 1e-8,            # xsf's documented large-k branch
    "zetac" => 1e-12,
    "trigamma" => 5e-14,
    "digammainv" => 1e-11,
    "gen_harmonic" => 1e-11,
    "stirling2" => 1e-5,        # xsf's Temme approximation at n = 220
)

# The double-double path is checked against mpmath. 5e-27 leaves three digits
# of margin over the ~1e-29 tolerances asserted by tenferro-sci-rs. Results
# below about 1e-292 have a subnormal low word, so there the relative accuracy
# falls towards `Float64`'s (see the `specials` crate documentation); such a
# result is checked at the `Float64` level instead.
const DD_TOL = 5e-27
const DD_SUBNORMAL_TOL = 1e-14
dd_tolerance(want) = abs(want) < big"1e-292" ? DD_SUBNORMAL_TOL : DD_TOL

call_f64(name, a, b) = b === nothing ? getfield(TenferroSci, Symbol(name))(a) :
                      getfield(TenferroSci, Symbol(name))(a, b)
call_dd(name, a, b) = b === nothing ? getfield(TenferroSci, Symbol(name))(Df64(a)) :
                     getfield(TenferroSci, Symbol(name))(Df64(a), Df64(b))
call_mf2(name, a, b) = b === nothing ? getfield(TenferroSci, Symbol(name))(Float64x2(a)) :
                      getfield(TenferroSci, Symbol(name))(Float64x2(a), Float64x2(b))

function relative_error(got::BigFloat, want::BigFloat)
    want == 0 && return iszero(got) ? BigFloat(0) : BigFloat(Inf)
    return abs(got - want) / abs(want)
end

@testset "Float64 against mpmath" begin
    for row in ROWS
        got = call_f64(row.name, row.a, row.b)
        want = row.value
        @test isfinite(got) || !isfinite(want)
        err = relative_error(BigFloat(got), want)
        @test err <= F64_TOL[row.name]
    end
end

@testset "Df64 against mpmath" begin
    for row in ROWS
        got = call_dd(row.name, row.a, row.b)
        want = row.value
        err = at_big_precision() do
            relative_error(big(got), want)
        end
        @test err <= dd_tolerance(want)
    end
end

@testset "special values" begin
    @test isnan(gamma(-1.0))
    @test gamma(0.0) == Inf
    @test gamma(-0.0) == -Inf
    @test lgamma(0.0) == Inf
    @test gammasgn(-0.5) == -1.0
    @test gammasgn(-1.5) == 1.0
    @test rgamma(-3.0) == 0.0
    @test zetac(1.0) == Inf
    @test zetac(Inf) == 0.0
    @test digamma(0.0) == -Inf
    @test digamma(-0.0) == Inf
    @test stirling2(5.0, 0.0) == 0.0
    @test stirling2(0.0, 0.0) == 1.0
    @test gen_harmonic(4.0, -1.0) == 10.0
    @test binom(10.0, 3.0) == 120.0

    @test isnan(gamma(Df64(-1.0)))
    @test isinf(gamma(Df64(0.0)))
    @test zetac(Df64(1.0)) == Df64(Inf)
    @test gammasgn(Df64(-0.5)) == Df64(-1.0)
    @test rgamma(Df64(-3.0)) == Df64(0.0)
end

@testset "Julia-facing argument types" begin
    @test gamma(5) == 24.0
    @test gamma(Float32(2.5)) == gamma(2.5)
    @test zeta(2, 1) == zeta(2.0, 1.0)
    @test stirling2(5, 3) == 25.0

    # BigFloat is not silently rounded to Float64.
    @test_throws MethodError gamma(big"2.5")
    @test big(Df64(big"2.5")) == 2.5

    # Mixed real / Df64 promotes to the double-double path.
    @test zeta(Df64(2.0), 1) == zeta(Df64(2.0), Df64(1.0))
    @test beta(Df64(2.0), 3.0) == beta(Df64(2.0), Df64(3.0))
    @test poch(3.0, Df64(4.0)) == poch(Df64(3.0), Df64(4.0))
end

@testset "Df64" begin
    d = Df64(2.5)
    @test d.hi == 2.5
    @test d.lo == 0.0
    @test big(d) == 2.5
    @test Float64(d) == 2.5
    @test isfinite(d) && !isinf(d) && !isnan(d) && !iszero(d)
    @test signbit(Df64(-1.0, 0.0))

    # A two-sum normalises an arbitrary finite pair.
    e = Df64(1.0, 1e-20)
    @test e.hi == 1.0
    @test e.lo == 1e-20

    # Non-finite values keep a zero low word.
    @test Df64(NaN).lo == 0.0
    @test Df64(Inf).lo == 0.0
    @test Df64(1.0, Inf) == Df64(Inf)
    @test Df64(floatmax(), floatmax()) == Df64(Inf)

    # 1/3 needs more than 53 bits; the double-double holds about 32 digits.
    third = Df64(1 // 3)
    @test abs(big(third) - big(1) / 3) / (big(1) / 3) < 1e-30
end

@testset "Float64x2 interface" begin
    for row in ROWS
        got = call_mf2(row.name, row.a, row.b)
        @test got isa Float64x2
        err = at_big_precision() do
            relative_error(BigFloat(got), row.value)
        end
        @test err <= dd_tolerance(row.value)
    end

    # The two-limb representations convert exactly.
    for (hi, lo) in ((1.0, 1e-20), (2.5, 0.0), (1e-300, 1e-320))
        d = Df64(hi, lo)
        @test Df64(Float64x2(d)) === d
    end
    # The Float64x2 methods convert to Df64 internally and back, so their
    # result is bit-identical to converting the Df64 result.
    let d = gamma(Df64(2.5)), mf = gamma(Float64x2(2.5))
        @test getfield(mf, :_limbs) == (d.hi, d.lo)
    end
    @test isnan(Df64(Float64x2(Df64(NaN))))

    # Mixed arguments promote to the double-double path.
    @test zeta(Float64x2(2), 1) isa Float64x2
    @test beta(Float64x2(2), Df64(3)) isa Float64x2
    @test zeta(Df64(2), Float64x2(1)) isa Float64x2

    # Special values survive the conversion.
    @test isinf(gamma(Float64x2(0)))
    @test isnan(gamma(Float64x2(-1)))
end
