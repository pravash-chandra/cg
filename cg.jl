
# define A and b
A = rand(5, 5)
b = zeros(5)

function cg(A, b, iter)
        x = zeros(length(b))
        r = b[:]

        p = r[:]

        for i = 1:iter
                # step length 
                alpha = r' * r / (p' * A * p)

                # approx solution
                x = x + alpha * p

                r_new = r - alpha * A * p
                # @show r_new
                beta = r_new' * r_new / (r' * r)
                p = r_new + beta * p
                # update r
                r = r_new
        end
        x
end

cg(A, b, 10)
