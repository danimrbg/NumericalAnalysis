function solve_integral(f::Function, a::Float64, b::Float64, number_of_subintervals::Int64)

    if a == b
        return 0
    end

    h = (b - a) / (2 * number_of_subintervals)
    s1 = 0
    s2 = 0

    for i = 1:2*number_of_subintervals-1
        if i % 2 == 0
            s1 += f(a + i * h)
        else
            s2 += f(a + i * h)
        end
    end

    return (h / 3) * (f(a) + 2 * s2 + 4 * s1 + f(b))
end
