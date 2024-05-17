using LinearAlgebra: I

function inverse_matrix(A::Matrix{Float64})
    n = size(A)[1]
    
    B = Matrix{Float64}(undef, n, n)
    for j = 1:n
        b = 1.0 * I[1:n,j]
        B[:,j] = solve_system_with_scaled_partial_pivoting(A, b)
    end
    return B
end