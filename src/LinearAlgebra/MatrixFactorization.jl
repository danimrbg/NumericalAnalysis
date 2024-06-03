using LinearAlgebra: I

function matrix_factorization(A::Matrix{Float64})
    if size(A)[1] != size(A)[2]
        println("A matriz precisa ser quadrada.")
        return
    end

    n = size(A)[1]

    if isnothing(solve_system_with_scaled_partial_pivoting(A, Vector{Float64}(undef, n)))
        return
    end

    n_digits = 8

    Aext = copy(A)
    P = 1.0 * Matrix(I,n,n)

    for i = 1:n-1
        p = i
        msg_erro = false

        for k = i:n
            if round(Aext[k, i], digits = n_digits) != 0
                p = k
                break
            end
        end

        if p != i
            v = copy(P[p, :])
            P[p, :] = P[i, :]
            P[i, :] = v
        end

        for j = i+1:n
            m = Aext[j, i] / Aext[i, i]
            Aext[j, :] = Aext[j, :] - m * Aext[i, :]
        end
    end

    if round(Aext[n, n], digits = n_digits) == 0
        println("A matriz é singular.")
        return
    end

    x = Vector{Float64}(undef, n)
    x[n] = Aext[n, n+1] / Aext[n, n]

    for i = n-1:-1:1
        sum = 0
        for j = i+1:n
            sum += Aext[i, j] * x[j]
        end
        x[i] = (Aext[i, n+1] - sum) / Aext[i, i]
    end

    return x
end

# teste