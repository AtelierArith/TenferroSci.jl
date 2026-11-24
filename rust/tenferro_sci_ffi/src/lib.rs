//! Julia-facing wrappers for the `specials` crate of tenferro-sci-rs.
//!
//! RustCall.jl scans this file for `#[julia]` items. Every function is
//! concrete because `#[julia]` cannot bind a generic Rust function: the
//! generic `specials` API is instantiated here for `f64` and for
//! `xprec::Df64`.
//!
//! * `f(...)` takes and returns `f64` and calls the xsf port.
//! * `f_dd(hi, lo, ..., out)` takes the double-double arguments as
//!   `hi + lo` pairs of `f64` and writes the double-double result into
//!   `out: &mut Dd`, calling the precision-adaptive path. A return slot is
//!   used because RustCall's FFI contract does not pass a `#[julia]` struct
//!   back by value.
//!
//! Values follow `specials`' contract: domain errors and overflow are
//! returned as NaN or ±∞, with no error flag.

// The `#[julia]` macro generates public `extern "C"` entry points that
// dereference their raw-pointer arguments. Those wrappers are generated code;
// the hand-written functions below take no raw pointers.
#![allow(clippy::not_unsafe_ptr_arg_deref)]

use rustcall_julia_macros::julia;
use xprec::Df64;

/// A double-double value as the pair `hi + lo`, matching `xprec::Df64`.
///
/// Julia unpacks the two fields into `TenferroSci.Df64`. The pair must be
/// normalised (`hi + lo == hi` in `f64`); `TenferroSci.Df64` guarantees that
/// for every value it constructs.
#[julia]
pub struct Dd {
    pub hi: f64,
    pub lo: f64,
}

#[julia]
impl Dd {
    /// Construct the double-double `hi + lo`.
    #[julia]
    pub fn new(hi: f64, lo: f64) -> Self {
        Dd { hi, lo }
    }
}

impl Dd {
    /// Write a double-double result through the `&mut Dd` return slot.
    fn set(&mut self, x: Df64) {
        self.hi = x.hi();
        self.lo = x.lo();
    }
}

/// Rebuild a `Df64` from the pair Julia passed, following `specials::Real`
/// (`from_dd`): non-finite values have no low word to keep.
fn df64(hi: f64, lo: f64) -> Df64 {
    if hi.is_finite() {
        Df64::new_full(hi, lo)
    } else {
        Df64::new(hi)
    }
}

// ---------------------------------------------------------------------------
// One-argument functions
// ---------------------------------------------------------------------------

#[julia]
pub fn gamma(x: f64) -> f64 {
    specials::gamma(x)
}

#[julia]
pub fn gamma_dd(hi: f64, lo: f64, out: &mut Dd) {
    out.set(specials::gamma(df64(hi, lo)));
}

#[julia]
pub fn lgamma(x: f64) -> f64 {
    specials::lgamma(x)
}

#[julia]
pub fn lgamma_dd(hi: f64, lo: f64, out: &mut Dd) {
    out.set(specials::lgamma(df64(hi, lo)));
}

#[julia]
pub fn gammasgn(x: f64) -> f64 {
    specials::gammasgn(x)
}

#[julia]
pub fn gammasgn_dd(hi: f64, lo: f64, out: &mut Dd) {
    out.set(specials::gammasgn(df64(hi, lo)));
}

#[julia]
pub fn rgamma(x: f64) -> f64 {
    specials::rgamma(x)
}

#[julia]
pub fn rgamma_dd(hi: f64, lo: f64, out: &mut Dd) {
    out.set(specials::rgamma(df64(hi, lo)));
}

#[julia]
pub fn digamma(x: f64) -> f64 {
    specials::digamma(x)
}

#[julia]
pub fn digamma_dd(hi: f64, lo: f64, out: &mut Dd) {
    out.set(specials::digamma(df64(hi, lo)));
}

#[julia]
pub fn zetac(x: f64) -> f64 {
    specials::zetac(x)
}

#[julia]
pub fn zetac_dd(hi: f64, lo: f64, out: &mut Dd) {
    out.set(specials::zetac(df64(hi, lo)));
}

#[julia]
pub fn trigamma(x: f64) -> f64 {
    specials::trigamma(x)
}

#[julia]
pub fn trigamma_dd(hi: f64, lo: f64, out: &mut Dd) {
    out.set(specials::trigamma(df64(hi, lo)));
}

#[julia]
pub fn digammainv(y: f64) -> f64 {
    specials::digammainv(y)
}

#[julia]
pub fn digammainv_dd(hi: f64, lo: f64, out: &mut Dd) {
    out.set(specials::digammainv(df64(hi, lo)));
}

#[julia]
pub fn lgam1p(x: f64) -> f64 {
    specials::lgam1p(x)
}

#[julia]
pub fn lgam1p_dd(hi: f64, lo: f64, out: &mut Dd) {
    out.set(specials::lgam1p(df64(hi, lo)));
}

// ---------------------------------------------------------------------------
// Two-argument functions
// ---------------------------------------------------------------------------

#[julia]
pub fn zeta(s: f64, q: f64) -> f64 {
    specials::zeta(s, q)
}

#[julia]
pub fn zeta_dd(s_hi: f64, s_lo: f64, q_hi: f64, q_lo: f64, out: &mut Dd) {
    out.set(specials::zeta(df64(s_hi, s_lo), df64(q_hi, q_lo)));
}

#[julia]
pub fn beta(a: f64, b: f64) -> f64 {
    specials::beta(a, b)
}

#[julia]
pub fn beta_dd(a_hi: f64, a_lo: f64, b_hi: f64, b_lo: f64, out: &mut Dd) {
    out.set(specials::beta(df64(a_hi, a_lo), df64(b_hi, b_lo)));
}

#[julia]
pub fn lbeta(a: f64, b: f64) -> f64 {
    specials::lbeta(a, b)
}

#[julia]
pub fn lbeta_dd(a_hi: f64, a_lo: f64, b_hi: f64, b_lo: f64, out: &mut Dd) {
    out.set(specials::lbeta(df64(a_hi, a_lo), df64(b_hi, b_lo)));
}

#[julia]
pub fn poch(a: f64, m: f64) -> f64 {
    specials::poch(a, m)
}

#[julia]
pub fn poch_dd(a_hi: f64, a_lo: f64, m_hi: f64, m_lo: f64, out: &mut Dd) {
    out.set(specials::poch(df64(a_hi, a_lo), df64(m_hi, m_lo)));
}

#[julia]
pub fn binom(n: f64, k: f64) -> f64 {
    specials::binom(n, k)
}

#[julia]
pub fn binom_dd(n_hi: f64, n_lo: f64, k_hi: f64, k_lo: f64, out: &mut Dd) {
    out.set(specials::binom(df64(n_hi, n_lo), df64(k_hi, k_lo)));
}

#[julia]
pub fn gen_harmonic(n: f64, a: f64) -> f64 {
    specials::gen_harmonic(n, a)
}

#[julia]
pub fn gen_harmonic_dd(n_hi: f64, n_lo: f64, a_hi: f64, a_lo: f64, out: &mut Dd) {
    out.set(specials::gen_harmonic(df64(n_hi, n_lo), df64(a_hi, a_lo)));
}

#[julia]
pub fn stirling2(n: f64, k: f64) -> f64 {
    specials::stirling2(n, k)
}

#[julia]
pub fn stirling2_dd(n_hi: f64, n_lo: f64, k_hi: f64, k_lo: f64, out: &mut Dd) {
    out.set(specials::stirling2(df64(n_hi, n_lo), df64(k_hi, k_lo)));
}
