module Shape

using LinearAlgebra

function write_bzr(surf, filename)
    open(filename, "w") do f
        n, m = size(surf)[1:2] .- 1
        println(f, "$n $m")
        for i in 0:n, j in 0:m
            p = surf[i+1,j+1,1:3]
            println(f, "$(p[1]) $(p[2]) $(p[3])")
        end
    end
end

function read_bzr(filename)
    read_numbers(f, numtype) = map(s -> parse(numtype, s), split(readline(f)))
    open(filename) do f
        n, m = read_numbers(f, Int)
        surf = Array{Float64}(undef, n + 1, m + 1, 3)
        for i in 0:n, j in 0:m
            surf[i+1,j+1,:] = read_numbers(f, Float64)
        end
        surf
    end
end

function bernstein_all(n, u)
    coeff = [1.0]
    for j in 1:n
        saved = 0.0
        for k in 1:j
            tmp = coeff[k]
            coeff[k] = saved + tmp * (1.0 - u)
            saved = tmp * u
        end
        push!(coeff, saved)
    end
    coeff
end

function eval(surf, uv)
    n, m = size(surf)[1:2] .- 1
    coeff_u = bernstein_all(n, uv[1])
    coeff_v = bernstein_all(m, uv[2])
    result = [0, 0, 0]
    for i in 0:n, j in 0:m
        result += surf[i+1,j+1,1:3] * coeff_u[i+1] * coeff_v[j+1]
    end
    result
end

function cpts_du(surf)
    n, m = size(surf)[1:2] .- 1
    result = Array{Float64}(undef, n, m + 1, 3)
    for i in 0:n-1, j in 0:m
        result[i+1,j+1,:] = n * (surf[i+2,j+1,:] - surf[i+1,j+1,:])
    end
    result
end

function cpts_dv(surf)
    n, m = size(surf)[1:2] .- 1
    result = Array{Float64}(undef, n + 1, m, 3)
    for i in 0:n, j in 0:m-1
        result[i+1,j+1,:] = m * (surf[i+1,j+2,:] - surf[i+1,j+1,:])
    end
    result
end

function rotate_controls(surf, angle)
    R = [cos(angle) -sin(angle) 0
         sin(angle)  cos(angle) 0
         0           0          1]
    n, m = size(surf)[1:2]
    result = Array{Float64}(undef, n, m, 3)
    for i in 1:n, j in 1:m
        result[i,j,:] = R * surf[i,j,:]
    end
    result
end

# Surface 1 and Surface 2 are the same surface with different parameterizations.
# Surface 3 is the same parameterization with Surface 2 but rotated in 3D space.
function test(newp = false)
    # triangle{1,2}[-new].bzr are two tensor product conversions
    # of a triangular Bezier patch, but with different parameterizations.
    # They evaluate to the same point at (u,v) = (0.5,0.5).
    local s
    if newp
        s = [read_bzr("triangle1-new.bzr"), read_bzr("triangle2-new.bzr")]
    else
        s = [read_bzr("triangle1.bzr"), read_bzr("triangle2.bzr")]
    end
    push!(s, rotate_controls(s[2], deg2rad(13)))
    uv = [0.5, 0.5]
    for i in 1:3
        x = eval(s[i], uv)
        cu, cv = cpts_du(s[i]), cpts_dv(s[i])
        xu, xv = eval(cu, uv), eval(cv, uv)
        n = normalize!(cross(xu, xv))
        xuu = eval(cpts_du(cu), uv)
        xuv = eval(cpts_dv(cu), uv)
        xvv = eval(cpts_dv(cv), uv)
        E, F, G = dot(xu, xu), dot(xu, xv), dot(xv, xv)
        L, M, N = dot(n, xuu), dot(n, xuv), dot(n, xvv)
        I = [E F; F G]
        II = [L M; M N]
        S = inv(I) * II
        J = inv(I) * [xu xv]'
        W = J' * II * J
        println("=== Surface $i")
        @show x
        @show n
        println("--- Shape operator (Weingarten map)")
        @show S
        @show det(S)
        @show tr(S) / 2
        eig = eigen(S)
        k1, k2 = eig.values
        @show (k1, k2)
        v = eig.vectors
        t1, t2 = v[:,1], v[:,2]
        @show (t1, t2)
        e1 = normalize!(v[1,1] * xu + v[2,1] * xv)
        e2 = normalize!(v[1,2] * xu + v[2,2] * xv)
        @show (e1, e2)
        println("--- Curvature tensor (Embedded Weingarten map)")
        @show W
        @show (tr(W)^2 - tr(W^2)) / 2
        @show tr(W) / 2
        eig = eigen(W)
        @show eig.values
        @show eig.vectors
    end
end

end
