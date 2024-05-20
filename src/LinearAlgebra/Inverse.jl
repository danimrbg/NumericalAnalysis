using LinearAlgebra: I

function inverse_matrix(A::Matrix{Float64})
    if size(A)[1] != size(A)[2]
        println("A matriz precisa ser quadrada.")
        return
    end

    n = size(A)[1]    
    B = Matrix{Float64}(undef, n, n)
    for j = 1:n
        b = 1.0 * I[1:n,j]
        if isnothing(solve_system_with_scaled_partial_pivoting(A, b))
            return
        end
        B[:,j] = solve_system_with_scaled_partial_pivoting(A, b)
    end
    return B
end