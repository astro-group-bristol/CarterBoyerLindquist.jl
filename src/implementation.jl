@with_kw struct CarterMethodBL{T} <: AbstractMetricParams{T}
    @deftype T
    M = 1.0
    a = 0.0
    E = 1.0
end

inner_radius(m::CarterMethodBL{T}) where {T} = m.M + √(m.M^2 - m.a^2)
constrain(::CarterMethodBL{T}, u, v; μ = 0.0) where {T} = v[1]

@with_kw struct CarterGeodesicPoint{T,P} <: AbstractGeodesicPoint{T}
    retcode::Symbol
    t::T
    u::AbstractVector{T}
    v::AbstractVector{T}
    p::P
end

function geodesic_point_type(m::CarterMethodBL{T}) where {T}
    p_type = typeof(make_parameters(T(0), T(0), 1))
    CarterGeodesicPoint{T,p_type}
end

make_parameters(L, Q, sign_θ) =
    (L = L, Q = Q, r = -1, θ = convert(Int, sign_θ), changes = T[0.0, 0.0])

function integrator_problem(
    m::CarterMethodBL{T},
    pos::StaticVector{S,T},
    vel::StaticVector{S,T},
    time_domain,
) where {S,T}
    L, Q = calc_lq(m, pos, vel)
    ODEProblem{false}(pos, time_domain, make_parameters(L, Q, vel[2])) do u, p, λ
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
    ODEProblem{true}(pos, time_domain, make_parameters(L, Q, vel[2])) do du, u, p, λ
        du .= carter_velocity(u, m.E, m.M, m.a, p)
    end
end

function metric_callback(m::CarterMethodBL{T}) where {T}
    (
        DiscreteCallback(ensure_domain(m), terminate!),
        DiscreteCallback(radial_negative_check(m), flip_radial_sign!),
        DiscreteCallback(angular_negative_check(m), flip_angular_sign!),
    )
end

function alpha_beta_to_vel(m::CarterMethodBL{T}, u, α, β) where {T}
    sinΦ, sinΨ = sinΦsinΨ(m.M, u[2], m.a, u[3], α, β)
    (β < 0.0 ? -1.0 : 1.0, sinΦ, sinΨ)
end

convert_velocity_type(u::StaticVector{S,T}, v) where {S,T} = convert(SVector{S,T}, v)
convert_velocity_type(u::AbstractVector{T}, v) where {T} = convert(typeof(u), collect(v))

function get_endpoint(
    m::CarterMethodBL{T},
    sol::SciMLBase.AbstractODESolution{T,N,S},
) where {T,N,S}
    us, ts, p = unpack_solution(sol)
    u = us[end]
    v = carter_velocity(u, m.E, m.M, m.a, p)
    t = ts[end]
    CarterGeodesicPoint(sol.retcode, t, u, convert_velocity_type(u, v), p)
end

export CarterMethodBL, CarterGeodesicPoint
