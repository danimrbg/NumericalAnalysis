function bisection(f::Function, a, b, τ::Float64 = 1e-15, N::Int64 = 10000)
    ai = a
    bi = b

    for i = 1:N
        p = (ai + bi) / 2
        fp = f(p)

        if fp ≈ 0
            atol = 1e-15 || (bi - ai) / 2 < τ
            return p
        end

        if f(ai) * fp > 0
            ai = p
        else
            bi = p
        end
    end

    return "O método falhou após " * string(N) * " interações."
end
