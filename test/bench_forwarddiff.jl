using ForwardDiff
using BenchmarkTools

const N = 100
const ITER = 1000

function build_inputs(n::Int)
    x = 0.01 .* collect(1:n)
    A = Array{Float64}(undef, n, n)
    @inbounds for j in 1:n
        for i in 1:n
            A[i, j] = 0.001 * sin(0.01 * (i + j)) + 0.001 * cos(0.02 * (i * j))
        end
    end
    return A, x
end

function dense_fun(x, A)
    @inbounds return sin.(x) .+ A * x
end

function run_bench(n::Int, iter::Int)
    A, x0 = build_inputs(n)
    f = let A = A
        x -> dense_fun(x, A)
    end

    ForwardDiff.jacobian(f, x0) # warm-up to avoid compilation time

    total = @belapsed begin
        for _ in 1:$iter
            ForwardDiff.jacobian($f, $x0)
        end
    end
    avg = total / iter

    println("n=$(n), iterations=$(iter), avg_time_s=$(avg)")
end

run_bench(N, ITER)
