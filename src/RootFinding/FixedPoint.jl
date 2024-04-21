function fixed_point(f, p0, τ::Float64 = 1e-15, N::Int64 = 10000)
    p = p0
    pi = p0

    for i = 1:N
        p = f(pi)

        if abs(p - pi) < τ
            return p
        end

        pi = p
    end
    
    return "O método falhou após " * string(N) * " interações."
end