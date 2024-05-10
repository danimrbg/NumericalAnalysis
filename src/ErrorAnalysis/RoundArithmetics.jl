function round_sum(x::Float64, y::Float64; digits::Int64)
    return round(round(x, digits = digits) + round(y, digits = digits), digits = digits)
end

function round_mul(x::Float64, y::Float64; digits::Int64)
    round(round(x, digits = digits) * round(y, digits = digits), digits = digits)
end

function trunc_sum(x::Float64, y::Float64; digits::Int64)
    return trunc(trunc(x, digits = digits) + trunc(y, digits = digits), digits = digits)
end

function trunc_mul(x::Float64, y::Float64; digits::Int64)
    trunc(trunc(x, digits = digits) * trunc(y, digits = digits), digits = digits)
end
