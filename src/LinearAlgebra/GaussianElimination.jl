function dimension_failure(A::Matrix{Float64}, b::Vector{Float64})
    if length(b) != size(A)[1]
        println("O tamanho da matriz e do vetor são incompatíveis.")
        return true
    end
    if size(A)[1] != size(A)[2]
        println("A matriz precisa ser quadrada.")
        return true
    end
    return false
end

function solve_system(A::Matrix{Float64}, b::Vector{Float64})
    dimension_failure(A, b) && return

    n_digits = 8

    Aext = [copy(A) copy(b)]
    n = length(b)

    for i = 1:n-1
        p = i
        msg_erro = false

        for k = i:n
            if round(Aext[k, i], digits = n_digits) != 0
                p = k
                break
            end

            if k == n
                msg_erro = true
            end
        end

        if msg_erro
            println("A matriz é singular.")
            return
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
