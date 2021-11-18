"""
Inputs:
    X: Historical demand values
    n: Number of periods to forecast forward
    a: Alpha smoothing value (0-1)
Outputs:
    F: Forecasted Demand Value
    M: Estimated demand for the fitted line
    S: Estimate of demand volume
    I: estiamte of intervals between demand

M forecasted demand values for the line fit. F is the last and therefore future forecasted value
q interval between last 2 periods with demand
S estimate of demand volume
I estimate of intervals between demand
"""

function CrostonSBA(X,n,alpha,disp_output=false)
    #needed if you want to plot the result
    
    q = 1
    a = alpha[1] #inputted as a vector/array and needs to be a number for calculations
    
    #initiate variables
    X_copy = copy(X)   #needs to happen so X is not modified
    S = zeros(size(X))
    I = zeros(size(X))
    
    #Initialize the procedure for the first observation
    if X_copy[1] > 0
        S[1] = X_copy[1]
        I[1] = 1
    else 
        # *** THIS IS NON STANDARD** read how to do this if first obs is zero
        X_copy[1] = 1
        S[1] = X_copy[1]
        I[1] = 1
    end
    
    if length(X) > 1
        for i in 2:length(X);
            if X_copy[i] == 0;
                S[i] = S[i-1];
                I[i] = I[i-1];        
                q = q + 1;
            else
                S[i]= a*X_copy[i].+(1-a)*S[i-1];
                I[i]= a*q.+((1-a)*I[i-1]);
                q = 1;
            end
        end
    end
    
    S = prepend!(S, NaN) # offset the vector so the time periods make logical sense.
    I = prepend!(I, NaN) # ie. now the s/i combo for period 2 is the forecast for period 2. Previously, the S/I combo for period 1 was the forecast for period 2.
    M = zeros(size(S)) #initiate M variable
    
    # Calculate Forecast
    for i in 1:length(S);
        M[i] = (1-a/2)*(S[i]/I[i]); # M is demand forecast over time for the associated time period. M25 is the forecasted value for period 25.
    end         # the last value of M is your new forecasted value

    F = M[end] # Forecasted Demand Value. (future)

    M = M[1:end-1] # This is your fitted line forecasted values. To be used for error calcs (in sample).
    
    if disp_output == true
        println("The forecasted Demand is: ", string(F), " for the next ", string(n), " periods.")
    end
    # Show plot of data
    #bar(X, xlabel="Time(days)", ylabel="Demand", title = "Demand Forecast (SBA)", label="Historical Demand")
    #plot!(M, label = "Fitted Line") #plot fitted line
    #plot!([length(M)+1, length(M)+1+60],[F, F], label="Forecasted") # Plot the forecast line
    
    return [F, M, S, I]

    # 
    # R Code that this is verified with
    # library(tsintermittent)
    # install.packages("readxl")
    # library("readxl")
    # my_data <- read_excel(file.choose())
    # x <- my_data[2]
    # crost(y,h=1,w=0.2,init="naive",type="sba",init.opt=FALSE, outplot = TRUE)  

end