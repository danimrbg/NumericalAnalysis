function newton(f, df, p0, τ::Float64 = 1e-15, N::Int64 = 10000)
    pi = p0
    p = p0

    for i = 1:N
        if df(pi) ≈ 0 atol = 1e-15
            return "Método falhou, pois a derivada anulou."
        end

        p = pi - f(pi) / df(pi)

        if abs(p-pi) < τ
            return p
        end

        pi = p
    end
    
    return "O método falhou após " * string(N) * " interações."
end