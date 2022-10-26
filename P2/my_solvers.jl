include("my_solvers_old.jl")




function conj_grad(A, b, guess, ϵ, iterMax)
    r = A*guess-b
    ρTwo = 0
    pOne = 0
    for i in 1:iterMax
        ρOne = transpose(r)*r
        if i == 1
            p = r
        else
            β = ρOne/ρTwo
            p = r+pOne*β    
        end
        q = A*p
        d = ρOne / (transpose(p)*q)
        guess = guess - p*d
        r = r - q*d
        ρTwo = ρOne
        pOne = p
        error = sqrt(transpose(A*guess-b)*(A*guess-b))/sqrt(transpose(guess)*guess)
        if error[1] <= ϵ
            print(i)
            print("\n")
            return guess
        end
    end
    return guess
end


