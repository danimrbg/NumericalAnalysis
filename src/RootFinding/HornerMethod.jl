function horner_method(a::Tuple, x0)
    n = length(a) - 1

    y = a[n+1]
    z = a[n+1]

    for j = n-1:-1:1
        y = x0 * y + a[j+1]
        z = x0 * z + y
    end

    y = x0 * y + a[1]

    return y, z
end
