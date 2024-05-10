module NumericalAnalysis

include("Includes.jl")

export round_sum,
    round_mul,
    trunc_sum,
    trunc_mul,
    absolute_error,
    relative_error,
    approximated_value_interval_given_relative_error,
    solve_system,
    bisection,
    fixed_point,
    newton,
    secant,
    delta_squared_aitken,
    steffensen,
    horner_method,
    solve_integral

end
