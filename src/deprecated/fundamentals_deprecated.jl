# Doesn't really work because indices change order in ITensors
function MPStoVector(mps::MPS)
    N = length(mps)
    A = [Array(mps[i], mps[i].inds...) for i=1:N]
    dims = size(A[1])
    A[1] = reshape(A[1], 1, dims...)
    A[1] = permutedims(A[1], [1,3,2])
    for i=2:N-1
        A[i] = permutedims(A[i], [3,1,2])
    end
    dims = size(A[N])
    A[N] = reshape(A[N], dims..., 1)
    A[N] = permutedims(A[N], [1,3,2])
    return A
end

