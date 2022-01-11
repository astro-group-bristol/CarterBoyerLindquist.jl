"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.10):

```math
T = E \\left( r^2 + a^2 \\right) - L * a.
```
"""
@inline T(E, L, r, a) = E * (r^2 + a^2) - L * a


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.10), for the special case of a null-geodesic ``\\mu = 0``:

```math
V_r = T^2 - \\Delta \\left[ (L - a E)^2 + Q \\right]
```
"""
@inline Vr(E, L, M, Q, r, a) = T(E, L, r, a)^2 - Δ(M, r, a) * ((L - a * E)^2 + Q)


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.10), for the special case of a null-geodesic ``\\mu = 0``:

```math
V_\\theta = 
    Q + \\cos^2 (\\theta) \\left[ a^2 E^2 - \\frac{L^2}{\\sin^2 (\\theta) } \\right].
```
"""
@inline Vθ(E, L, Q, a, θ) = Q + cos(θ)^2 * ((a * E)^2 - (L / sin(θ))^2)


"""
    $(TYPEDSIGNATURES)

The ``t`` compontent of the equation of motion for a photon around a black hole, multiplied 
by ``\\Sigma``.

From Bardeen et al. (1972) eq. (2.9d):

```math
\\Sigma \\frac{\\text{d}t}{\\text{d}\\lambda} =
- a \\left( a E \\sin^2 \\theta - L \\right)
+ \\frac{\\left( r^2 + a^2 \\right) T}{\\Delta}.
```
"""
@inline function Σδt_δλ(E, L, M, r, a, θ)
    -a * (a * E * sin(θ)^2 - L) + (r^2 + a^2) * T(E, L, r, a) / Δ(M, r, a)
end


"""
    $(TYPEDSIGNATURES)

The ``r`` compontent of the equation of motion for a photon around a black hole, multiplied 
by ``\\Sigma``.

Modified from Bardeen et al. (1972) eq. (2.9a):

```math
\\Sigma \\frac{\\text{d}r}{\\text{d}\\lambda} =
\\pm \\sqrt{\\lvert V_r \\rvert},
```

where, for implementation reason, the sign is always positive. Instead, the sign is applied 
in [`δ`](@ref).
"""
@inline function Σδr_δλ(E, L, M, Q, r, a)
    V = Vr(E, L, M, Q, r, a)
    √(V * sign(V))
end


"""
    $(TYPEDSIGNATURES)

The ``\\theta`` compontent of the equation of motion for a photon around a black hole, 
multiplied by ``\\Sigma``.

Modified from Bardeen et al. (1972) eq. (2.9b):

```math
\\Sigma \\frac{\\text{d}\\theta}{\\text{d}\\lambda} =
\\pm \\sqrt{\\lvert V_\\theta \\rvert},
```

where, for implementation reason, the sign is always positive. Instead, the sign is applied 
in [`δ`](@ref).
"""
@inline function Σδθ_δλ(E, L, Q, a, θ)
    V = Vθ(E, L, Q, a, θ)
    √(V * sign(V))
end


"""
    $(TYPEDSIGNATURES)

The ``\\phi`` compontent of the equation of motion for a photon around a black hole, 
multiplied by ``\\Sigma``.

From Bardeen et al. (1972) eq. (2.9c):

```math
\\Sigma \\frac{\\text{d}\\phi}{\\text{d}\\lambda} =
- \\frac{L}{\\sin^2 \\theta} - aE
+ \\frac{aT}{\\Delta}.
```
"""
@inline function Σδϕ_δλ(E, L, M, r, a, θ)
    (L / sin(θ)^2) - (a * E) + a * T(E, L, r, a) / Δ(M, r, a)
end


function carter_velocity(u, E, M, a, p)
    let r = u[2], θ = u[3], L = p.L, Q = p.Q
        Σ₀ = Σ(r, a, θ)
        @show (
            Σδt_δλ(E, L, M, r, a, θ) / Σ₀,
            p.r * Σδr_δλ(E, L, M, Q, r, a) / Σ₀,
            p.θ * Σδθ_δλ(E, L, Q, a, θ) / Σ₀,
            Σδϕ_δλ(E, L, M, r, a, θ) / Σ₀
        )
    end
end