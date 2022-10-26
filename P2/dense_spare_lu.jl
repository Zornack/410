using LinearAlgebra
using SparseArrays
function denseMatrix()
    N = 1000
    A = zeros(N,N)
    b = ones(N)
    for i in 1:N
        A[i,i] = -2
        if i == 1
            A[i,i+1] = 1
        elseif i == N
            A[i,i-1] = 1
        else
            A[i,i+1] = 1
            A[i,i-1] = 1
        end
    end
    t1 = @elapsed F = lu(A)
    t1 = @elapsed F = lu(A)
    t2 = @elapsed x = F\b
    t2 = @elapsed x = F\b
    return t1,t2
end



function sparseMatrix()
    N = 1000
    A = spzeros(N,N)
    b = ones(N)
    for i in 1:N
        A[i,i] = -2
        if i == 1
            A[i,i+1] = 1
        elseif i == N
            A[i,i-1] = 1
        else
            A[i,i+1] = 1
            A[i,i-1] = 1
        end
    end
    t1 = @elapsed F = lu(A)
    t1 = @elapsed F = lu(A)
    t2 = @elapsed x = F\b
    t2 = @elapsed x = F\b
    return t1,t2
end

N = 1000
Ad = zeros(N,N)
bd = ones(N)
for i in 1:N
    Ad[i,i] = -2
    if i == 1
        Ad[i,i+1] = 1
    elseif i == N
        Ad[i,i-1] = 1
    else
        Ad[i,i+1] = 1
        Ad[i,i-1] = 1
    end
end
t1d = @elapsed Fd = lu(Ad)
t1d2 = @elapsed Fd2 = lu(Ad)
t2d = @elapsed xd = Fd\bd
t2d2 = @elapsed xd2 = Fd\bd

As = spzeros(N,N)
bs = ones(N)
for i in 1:N
    As[i,i] = -2
    if i == 1
        As[i,i+1] = 1
    elseif i == N
        As[i,i-1] = 1
    else
        As[i,i+1] = 1
        As[i,i-1] = 1
    end
end
t1s = @elapsed Fs = lu(As)
t1s2 = @elapsed Fs2 = lu(As)
t2s = @elapsed xs = Fs\bs
t2s2= @elapsed xs2 = Fs\bs