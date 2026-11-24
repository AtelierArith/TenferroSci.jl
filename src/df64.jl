"""
    Df64(hi, lo)
    Df64(x)

A double-double number, the pair `hi + lo` of two `Float64`s, matching the
`xprec::Df64` type used by the Rust side (about 32 significant digits).

`Df64(x)` converts a real number. A `Float64` is exact; other reals go through
`BigFloat`, so a value needing more than the 106 bits of a double-double is
rounded. `Df64(hi, lo)` is normalised with a two-sum, so any pair of finite
`Float64`s is accepted. A non-finite `hi` keeps a zero low word, as in
`xprec`.

`Df64` is deliberately not a subtype of `Real`: it supports the special
functions of this package, conversion with [`Base.big`](@ref) and
[`Float64`](@ref), and inspection (`isnan`, `isinf`, `isfinite`, `iszero`,
`signbit`), but no arithmetic. Do arithmetic in `BigFloat` or `Float64` and
convert the result.

```julia
julia> d = Df64(2.5);

julia> big(d)
2.5

julia> Float64(d)
2.5

julia> Df64(1.0, 1e-20).lo
1.0e-20
```
"""
struct Df64
    hi::Float64
    lo::Float64

    function Df64(hi::Float64, lo::Float64)
        if isfinite(hi) && isfinite(lo)
            # Knuth's two-sum: s + e == hi + lo exactly, and s + e == s.
            s = hi + lo
            if isfinite(s)
                bp = s - hi
                e = (hi - (s - bp)) + (lo - bp)
                return new(s, e)
            end
            return new(s, 0.0)   # the pair overflows the f64 range
        elseif isfinite(hi)
            return new(hi + lo, 0.0)  # only `lo` is non-finite
        else
            return new(hi, 0.0)
        end
    end
end

Df64(x::Df64) = x
Df64(x::Float64) = Df64(x, 0.0)
Df64(x::Real) = Df64(BigFloat(x))

function Df64(x::BigFloat)
    hi = Float64(x)
    isfinite(hi) || return Df64(hi, 0.0)
    lo = Float64(x - hi)
    Df64(hi, lo)
end

Base.big(d::Df64) = BigFloat(d.hi) + BigFloat(d.lo)
Base.Float64(d::Df64) = d.hi
Base.float(d::Df64) = d

Base.isfinite(d::Df64) = isfinite(d.hi)
Base.isinf(d::Df64) = isinf(d.hi)
Base.isnan(d::Df64) = isnan(d.hi)
Base.iszero(d::Df64) = iszero(d.hi) && iszero(d.lo)
Base.signbit(d::Df64) = signbit(d.hi)

Base.zero(::Type{Df64}) = Df64(0.0, 0.0)
Base.one(::Type{Df64}) = Df64(1.0, 0.0)

function Base.show(io::IO, d::Df64)
    print(io, "Df64(", d.hi, ", ", d.lo, ")")
end

"""Unpack the Rust-owned `Bindings.Dd` handle into an immutable `Df64`."""
@inline unpack_dd(d) = Df64(d.hi, d.lo)
