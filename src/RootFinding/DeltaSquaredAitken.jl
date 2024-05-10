function delta_squared_aitken(f::Function, p0, τ::Float64 = 1e-15, N::Int64 = 10000)
    p1 = f(p0)
    p2 = f(p1)
    p = p2

    for i = 1:N
        if abs((p2 - p1) - (p1 - p0)) < 1e-16
            print("A sequência admite denominador muito próximo de zero.\n")
            return p
        end

        p = p0 - ((p1 - p0)^2 / ((p2 - p1) - (p1 - p0)))

        if abs(p - p2) < τ
            return p
        end

        p0 = p1
        p1 = p2
        p2 = f(p1)
    end

    return "O método falhou após " * string(N) * " interações."
end
