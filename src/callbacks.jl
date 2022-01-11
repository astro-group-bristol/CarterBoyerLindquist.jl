function ensure_domain(m::AbstractMetricParams{T}) where {T}
    min_r = inner_radius(m)
    (u, λ, integrator) -> u[2] ≤ min_r * 1.1 || u[2] > 1200.0
end

function flip_radial_sign!(integrator)
    p = integrator.p
    integrator.p = @set p.r = -p.r
end

function flip_angular_sign!(integrator)
    p = integrator.p
    integrator.p = @set p.θ = -p.θ
end

function radial_negative_check(m::CarterMethodBL{T}) where {T}
    (u, λ, integrator) -> begin
        L, Q, _, _ = integrator.p
        Vr(m.E, L, m.M, Q, u[2], m.a) < 0
    end
end

function angular_negative_check(m::CarterMethodBL{T}) where {T}
    (u, λ, integrator) -> begin
        L, Q, _, _ = integrator.p
        Vθ(m.E, L, Q, m.a, u[3]) < 0
    end
end