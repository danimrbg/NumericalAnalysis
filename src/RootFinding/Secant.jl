function secant(f::Function, p0, p1, τ::Float64 = 1e-15, N::Int64 = 10000)
    if abs(f(p0) - f(p1)) < 1e-15
        return "As imagens dos pontos dados estão próximas."
    end

    p = p1
    pi = p1
    qi = p0
    fpi = f(pi)
    fqi = f(qi)

    for i = 1:N
        if abs(fqi - fpi) < 1e-15
            return "Erro, pois a imagem de p na " *
                   string(i) *
                   "-ésima iteração está muito próxima da imagem de p na " *
                   string(i - 1) *
                   "-ésima iteração."
        end

        p = pi - (fpi * (pi - qi)) / (fpi - fqi)

        if abs(pi - qi) < τ
            return p
        end

        qi = pi
        pi = p
        fpi = f(pi)
        fqi = f(qi)
    end

    return "O método falhou após " * string(N) * " interações."
end