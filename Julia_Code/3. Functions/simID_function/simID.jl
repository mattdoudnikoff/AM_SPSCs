using Distributions
function simID(n, obs, adi, cv2, level)
    # n: number of demand streams to create
    # obs: number of time periods to simulate
    # adi: average demand interval
    # cv2: coefficient of variance squared
    # level: average of nonzero demand periods

    #d: 'n' x 'obs' matrix with the specified demand streams
    
    ########### ADI cannot be less than one. and cv2 cannot be < 0 
    #Error checking on ADI and CV2
    if adi < 1
        error("ERROR: ADI less than one. Pick value greater than or equal to 1")
    end
    if cv2 < 0
        error("ERROR: CV2 less than zero. Pick value greater than zero")
    end
    
    m = level - 1
    if (cv2 != 0)
        p = (m / ( cv2 * ( (m+1)^2) ) ); # calculates the p for negative binomial function
        r = m * p / (1 - p); # calculates the r for the negative binomial funcion
        if r<= 0
            error("NaNs Produced: Choose different cv2 and level combo")
        end
        d = rand(Binomial(1, 1/adi), n, obs) .* (rand(NegativeBinomial(r, p), n, obs) .+ 1)
        else 
        d = rand(Binomial(1, 1/adi), n, obs) .* (ones(5,6) * 11)
    end
    return d
end