function secant(f, p0, p1, τ::Float64 = 1e-15, N::Int64 = 10000)
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
            return "Erro, pois a imagem de p na " * string(i) * "-ésima iteração está muito próxima da imagem de p na " * string(i-1) * "-ésima iteração."
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

"""
INPUT: f, p0, p1, τ (tolerˆancia) e N (n´umero m´aximo de itera¸c˜oes);
OUTPUT: Aproxima¸c˜ao para p ou mensagem de erro;
STEP 1: if f(p0) = f(p1) then
    OUTPUT(“f(p0) e f(p1) est˜ao pr´oximos.”);
    STOP;
STEP 2: p := p1, pi := p1, qi := p0, fpi := f(pi), fqi := f(qi);
STEP 3: for i = 1, · · · , N (do STEPS 4-7):
    STEP 4: if fqi = fpi then
        OUTPUT(“Erro, pois f(pi) muito pr´oximo de f(pi)”);
        STOP;
    STEP 5: p := pi − fpi(pi − qi)fpi − fqi ;
    STEP 6: if |pi − qi| < τ then
        OUTPUT(p);
        STOP;
    STEP 7: qi := pi, pi := p; fpi := f(pi), fqi := f(qi)
"""