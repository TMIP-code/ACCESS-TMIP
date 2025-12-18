
# Ns = rand(5:50, 10)
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using CUDA
using LinearAlgebra
using SparseArrays

# Dummy block sizes
Ns = rand(3:5, 3)
N = sum(Ns)

# Dummy big Tridiagonal matrix
A = Tridiagonal(blockdiag([sparse(Tridiagonal(ones(N - 1), [1; 2 * ones(N - 2); 1], 1.1ones(N - 1))) + I for N in Ns]...))

# Block indices
blk_last_idx = cumsum(Ns);
blk_first_idx = [1; blk_last_idx[1:end-1] .+ 1];
blk_idx = [collect(first:last) for (first, last) in zip(blk_first_idx, blk_last_idx)]
offdiag_idx = [collect(first:(last - 1)) for (first, last) in zip(blk_first_idx, blk_last_idx)]

# Display blocks
for idx in blk_idx
    display(Tridiagonal(A[idx, idx]))
end

# Should I Preallocate big vectors on the GPU?
b = ones(N)
b_gpu = CuVector(b)


# Should I preallocate big vectors of the diagonals on the GPU?
dl_big = A.dl |> Vector |> CuVector
d_big = A.d |> Vector |> CuVector
du_big = A.du |> Vector |> CuVector

# Should I preallocate vectors of the blocks' first/last indices on the GPU
blk_last_idx_gpu = CuVector(blk_last_idx)
blk_first_idx_gpu = CuVector(blk_first_idx)

# Function to Solve each block on the GPU
function tridiagonalsolve!(dl_big, d_big, du_big, b_gpu, blk_first_idx_gpu, blk_last_idx_gpu)
    for (first, last) in zip(blk_first_idx_gpu, blk_last_idx_gpu)
        d = view(d_big, first:last)
        B = view(b_gpu, first:last)
        dl = view(dl_big, first:(last - 1))
        du = view(du_big, first:(last - 1))
        CUSPARSE.gtsv2!(dl, d, du, B)
    end
    return nothing
end

# # Solve on the GPU
x = copy(b)
x_gpu = CuVector(x)
@cuda tridiagonalsolve!(dl_big, d_big, du_big, x_gpu, blk_first_idx_gpu, blk_last_idx_gpu)

# Check correctness
x_cpu = A \ b
[Vector(x_gpu) x_cpu]


