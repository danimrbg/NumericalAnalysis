@testset "Test GaussianEliminationScaledPartialPivoting solve_system_with_scaled_partial_pivoting" begin

    atol = 1e-10

    A = [1.0 -1 1 1; 2 -1 -1 0; 1 -2 0 1; 5 0 1 -1]
    b = [4.0, -3, 1, 5]
    x = [0, -1, 4, -1]
    x_calc = solve_system_with_scaled_partial_pivoting(A, b)
    @test norm(x - x_calc, Inf) ≈ 0 atol = atol

    A = [1.0 -1 1 1; 2 -1 -1 0; 1 -2 0 1; 5 0 1 -1]
    b = [4.0, -3, 1, 4]
    x = [0, -2 / 3, 11 / 3, -1 / 3]
    x_calc = solve_system_with_scaled_partial_pivoting(A, b)
    @test norm(x - x_calc, Inf) ≈ 0 atol = atol

    A = [0.003 59.14; 5.291 -6.130]
    b = [59.17; 46.78]
    x = [10.0; 1.0]
    x_calc = solve_system_with_scaled_partial_pivoting(A, b)
    @test norm(x - x_calc, Inf) ≈ 0 atol = atol
end
