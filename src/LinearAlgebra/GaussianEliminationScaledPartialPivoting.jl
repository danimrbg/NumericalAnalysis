function solve_system_with_scaled_partial_pivoting(A::Matrix{Float64}, b::Vector{Float64})
    dimension_failure(A, b) && return

    n_digits = 8
    Aext = [copy(A) copy(b)]
    n = length(b)

    for i = 1:n-1
        p = i
        sk = Vector{Float64}(undef, n - i + 1)

        msg_erro = false

        for k = i:n
            sk[k-i+1] = 0

            for j = i:n
                if abs(Aext[k, j]) > sk[k-i+1]
                    sk[k-i+1] = abs(Aext[k, j])
                end
            end

            if round(sk[k-i+1], digits = n_digits) == 0
                msg_erro = true
            end
        end

        if msg_erro
            println("Sistema não admite única solução.")
            return
        end

        for k = i:n
            if abs(Aext[p, i]) / sk[p-i+1] < abs(Aext[k, i]) / sk[k-i+1]
                p = k
            end
        end

        if p != i
            v = copy(Aext[p, :])
            Aext[p, :] = Aext[i, :]
            Aext[i, :] = v
        end

        for j = i+1:n
            m = Aext[j, i] / Aext[i, i]
            Aext[j, :] = Aext[j, :] - m * Aext[i, :]
        end
    end

    if round(Aext[n, n], digits = n_digits) == 0
        println("Sistema não admite única solução.")
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
