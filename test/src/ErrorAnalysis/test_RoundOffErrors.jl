@testset "Test RoundOffErrors absolute_error" begin

    atol = 1e-15

    #Exercise 1_3a.
    value = Float64(π)
    approximated_value = 22 / 7

    @test absolute_error(value, approximated_value) ≈ 0.0012644892673496777 atol = atol

    #Exercise 1_3b.
    value = Float64(π)
    approximated_value = 3.1416

    @test absolute_error(value, approximated_value) ≈ 7.346410206832132e-6 atol = atol
end

@testset "Test RoundOffErrors relative_error" begin

    atol = 1e-15

    function stirling_approximation(n::Int64)::Float64
        sqrt(2 * pi * n) * (n / exp(1))^n
    end

    #Exercise 1_4a. 
    value = factorial(7) |> Float64
    approximated_value = stirling_approximation(7)

    @test relative_error(value, approximated_value) ≈ 0.011826223886416323 atol = atol

    #Exercise 1_4b. 
    value = factorial(20) |> Float64
    approximated_value = stirling_approximation(20)

    @test relative_error(value, approximated_value) ≈ 0.00415765262287991 atol = atol
end

@testset "Test RoundOffErrors approximated_value_interval_given_relative_error" begin

    atol = 1e-8

    println("************************************************")
    println("Remark: next test performs 2 tests per exercise!")
    println("************************************************")

    #Exercise 1_5a 
    min_approx, max_approx = approximated_value_interval_given_relative_error(10.0, 1e-3)

    @test min_approx ≈ 9.99 atol = atol
    @test max_approx ≈ 10.01 atol = atol

    #Exercise 1_5b
    min_approx, max_approx = approximated_value_interval_given_relative_error(100.0, 1e-3)

    @test min_approx ≈ 99.9 atol = atol
    @test max_approx ≈ 100.1 atol = atol

    #Exercise 1_5c 
    min_approx, max_approx = approximated_value_interval_given_relative_error(-50.0, 1e-3)

    @test min_approx ≈ -50.05 atol = atol
    @test max_approx ≈ -49.95 atol = atol

    #Exercise 1_5d
    min_approx, max_approx = approximated_value_interval_given_relative_error(10.0, 1e-4)

    @test min_approx ≈ 9.999 atol = atol
    @test max_approx ≈ 10.001 atol = atol

    #Exercise 1_5e
    min_approx, max_approx = approximated_value_interval_given_relative_error(100.0, 1e-4)

    @test min_approx ≈ 99.99 atol = atol
    @test max_approx ≈ 100.01 atol = atol

    #Exercise 1_5f
    min_approx, max_approx = approximated_value_interval_given_relative_error(-50.0, 1e-4)

    @test min_approx ≈ -50.005 atol = atol
    @test max_approx ≈ -49.995 atol = atol
end
