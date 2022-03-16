module CarterBoyerLindquist

using Accessors
using Parameters
using DocStringExtensions
using StaticArrays
using DifferentialEquations

import GeodesicBase:
    AbstractMetricParams,
    inner_radius,
    AbstractGeodesicPoint,
    get_endpoint,
    geodesic_point_type,
    unpack_solution,
    SciMLBase

import GeodesicTracer:
    DiscreteCallback,
    terminate!,
    integrator_problem,
    metric_callback,
    create_callback_set,
    constrain,
    alpha_beta_to_vel

include("coordinate-functions.jl")
include("motion-functions.jl")
include("implementation.jl")
include("motion-constants.jl")
include("callbacks.jl")

end # module
