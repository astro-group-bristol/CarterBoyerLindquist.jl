function carter_velocity(u, E, M, a, p)
    let r = u[2], θ = u[3], L = p.L, Q = p.Q
        Σ₀ = Σ(r, a, θ)
        (
            Σδt_δλ(E, L, M, r, a, θ) / Σ₀,
            p.r * Σδr_δλ(E, L, M, Q, r, a) / Σ₀,
            p.θ * Σδθ_δλ(E, L, Q, a, θ) / Σ₀,
            Σδϕ_δλ(E, L, M, r, a, θ) / Σ₀,
        )
    end
end

calc_lq(m::CarterMethodBL{T}, pos, vel) where {T} =
    LQ(m.M, pos[2], m.a, pos[3], vel[3], vel[4])


@with_kw struct CarterMethodBL{T} <: AbstractMetricParams{T}
    @deftype T
    M = 1.0
    a = 0.0
    E = 1.0
end

inner_radius(m::CarterMethodBL{T}) where {T} = m.M + √(m.M^2 - m.a^2)
constrain(::CarterMethodBL{T}, u, v; μ = 0.0) where {T} = v[1]


function integrator_problem(
    m::CarterMethodBL{T},
    pos::StaticVector{S,T},
    vel::StaticVector{S,T},
    time_domain,
) where {S,T}
    L, Q = calc_lq(m, pos, vel)
    ODEProblem{false}(pos, time_domain, make_parameters(L, Q, vel[2], T)) do u, p, λ
        SVector(carter_velocity(u, m.E, m.M, m.a, p)...)
    end
end

function integrator_problem(
    m::CarterMethodBL{T},
    pos::AbstractVector{T},
    vel::AbstractVector{T},
    time_domain,
) where {T}
    L, Q = calc_lq(m, pos, vel)
    ODEProblem{true}(pos, time_domain, make_parameters(L, Q, vel[2], T)) do du, u, p, λ
        du .= carter_velocity(u, m.E, m.M, m.a, p)
    end
end

export CarterMethodBL