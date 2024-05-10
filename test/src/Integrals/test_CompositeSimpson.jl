@testset "Test CompositeSimpson solve_integral" begin

    tol = 1.0e-3

    f(x, p) = sin(x * p)
    p = 1.7
    g(x) = sin(x * p)
    a = -2.0
    b = 5.0
    domain = (a, b)
    prob = IntegralProblem(f, domain, p)
    sol = solve(prob, QuadGKJL())

    # Função trigonométrica - caso 1
    @test solve_integral(g, a, b, 10_000) ≈ sol.u atol = tol

    # Função trigonométrica - caso 2
    @test solve_integral(g, b, a, 10_000) ≈ -sol.u atol = tol

    # Função trigonométrica - caso 3
    @test solve_integral(g, a, a, 1) == 0

    # Função polinomial
    tol = 1.0e-2

    f(x, p) = 2 * p * x^2 - 6
    p = 1.0
    g(x) = 2 * p * x^2 - 6
    a = 5.0
    b = 8.0
    domain = (a, b)
    prob = IntegralProblem(f, domain, p)
    sol = solve(prob, QuadGKJL())

    @test solve_integral(g, a, b, 10_000) ≈ sol.u atol = tol

    # Função exponencial

    tol = 1.0e-3

    f(x, p) = exp(x * p)
    p = 3.0
    g(x) = exp(x * p)
    a = -3.0
    b = 1.0
    domain = (a, b)
    prob = IntegralProblem(f, domain, p)
    sol = solve(prob, QuadGKJL())

    @test solve_integral(g, a, b, 100_000) ≈ sol.u atol = tol

    # Função logarítmica

    tol = 1.0e-3

    f(x, p) = log(exp(1), x * p)
    p = 5.0
    g(x) = log(exp(1), x * p)
    a = 2.0
    b = 9.0
    domain = (a, b)
    prob = IntegralProblem(f, domain, p)
    sol = solve(prob, QuadGKJL())

    @test solve_integral(g, a, b, 10_000) ≈ sol.u atol = tol
end
