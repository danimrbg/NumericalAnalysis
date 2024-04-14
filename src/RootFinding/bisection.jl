function bisection(f, a, b, τ::Float64 = 1e-5, N::Int64 = 100)
    ai = a
    bi = b
    for i in 1:N
        p = (ai + bi) / 2
        fp = f(p)
        if fp == 0 || (bi - ai) / 2 < τ
            return p
        end
        if f(ai) * fp > 0
            ai = p
        else 
            bi = p
        end
    end
    return "O método falhou."
end