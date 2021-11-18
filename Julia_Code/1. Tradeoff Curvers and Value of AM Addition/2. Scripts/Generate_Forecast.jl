###needed referred files
#This work is done to automate transfer runs between computers and OSs (only Mac and Windows are supported)
file_path = pwd()
if file_path[1] == '/'
    inipath = file_path[1:findlast("Julia_Code",file_path)[end]]*"/"
    include(inipath * "3. Functions/simID_function/simID.jl")
    include(inipath * "3. Functions/simID_function/idClass.jl")
    include(inipath * "3. Functions/Forecasting_Functions/CrostonSBA.jl")
    include(inipath * "3. Functions/Forecasting_Functions/Croston.jl")
    include(inipath * "3. Functions/Forecasting_Functions/MarkovBootstrap.jl")
else
    inipath = file_path[1:findlast("Julia_Code",file_path)[end]]*"\\"
    include(inipath * "3. Functions\\simID_function\\simID.jl")
    include(inipath * "3. Functions\\simID_function\\idClass.jl")
    include(inipath * "3. Functions\\Forecasting_Functions\\CrostonSBA.jl")
    include(inipath * "3. Functions\\Forecasting_Functions\\Croston.jl")
    include(inipath * "3. Functions\\Forecasting_Functions\\MarkovBootstrap.jl")
end

#Packages needed
using Optim
using Plots
using StatsBase
using Distributions
using Statistics
using Random

"""
    Generate historical demand per product for each distribution center. 
Inputs:
    adi: 
    CV2: 
    lev: 
    G: Number of products
    ncust: Number of customers
    obs: Number of time periods of historical data to generate
Outputs: 
    (ncust x obs x G) array of historical demand
"""
function gen_historical(adi, cv2, lev, G, nsup, obs=180)
    if G < 1 
        error("ERROR: G is less than 1, pick value of 1 or more")
    end
    if nsup < 0
        error("ERROR: ncust is less than 1, pick value of 1 or more")
    end
    
    historical = Array{Int64}(undef,nsup,obs,G)
    
    for g in 1:G #iterate for number of products
        for ns in 1:nsup #iterate for number of suppliers
            historical[ns,:,g] = simID(1,obs,adi,cv2,lev)
        end
    end
    return historical
end

"""
Generate demand realizations based off original input parameters (to be used in EVAL)
Inputs:
    nreal: number of realizations: 1st dimension
    adi:
    cv2:
    lev:
    T: Time horizon to be evaluated: dimension
    G: Number of products: 4th dimension
    nsup: Number of Suppliers: 3rd dimension
    dem_level: (OPTIONAL) When varying ADI and CV2, this is a nsup x G array that will keep the randomized data within 5% of this demand level for prouct and DC
Outputs:
""" ###MIGHT NEED TO EDIT THIS SLIGHTLY ESPECIALLY IF YOU'RE GOING TO VARY ADI AND CV2 @SEAN
function gen_dem_realization(nreal, adi, cv2, lev, T, G, nsup, dem_level=0)
    obs = T
    sup_dem = zeros(nreal,obs,nsup,G)
    
    ####generation procedure
    #Generate demand realization and determine the total qty of demand generated
    #Determine percent difference between baseline qty and generated qty
    #if difference is greater than desired, try again
    #if difference is within demand tolerance, accept demand realization.
    base_qty = zeros(nsup,G)
    for g = 1:G
        for ns = 1:nsup
            if dem_level == 0 #### this might need to be changed later
                #setting baseline to keep total qty of demand generated to within 5%
                baseline_demand = floor.(simID(20,T,adi,cv2,lev))  #20 realizations to ensure total demand is roughly in center of normal distribution of randomness
                base_qty[ns,g] = mean(sum(baseline_demand,dims=2)) #average total demand
            else
                base_qty[ns,g] = dem_level[ns,g]
            end
            for real = 1:nreal
                log = 1
                num_trys = 1
                raw_stream = simID(1,T,adi,cv2,lev) #initialize variable (would be forgotten otherwise because Julia's operating system)
                while log >= 0.0025 #5% gap between demand realizations (from base_qty)
                    num_trys = num_trys + 1
                    raw_stream = simID(1,T,adi,cv2,lev) #demand stream (a single potential realization)
                    raw_tot_qty = sum(raw_stream) # total qty generated
                    log = ((base_qty[ns,g] - raw_tot_qty)/base_qty[ns,g])^2 #difference between raw generated stream and the baseline dem qty
                    if num_trys == 2000 #if it's a weird draw and cannot get sample within tolerance
                        lev = lev + 1 #adjust level
                        num_trys = 1 #reset
                    end
                end
                sup_dem[real,:,ns,g] = raw_stream
            end
        end
    end
    
    return sup_dem, base_qty
end

"""
Function that first takes historical demand then classifies it and then gives a forecast
Inputs:
    hist_dem: historical demand that should already be generated
    G: number of Customers
    T: Time horizon interested in predicting
Outputs:
    A forecast of future demand
"""
function forecast_dem(hist_dem,G,T,ncust,DCrte)
    #loop through each product type of all customers
    nsup = size(hist_dem)[1]
    forecast = Array{Int64}(undef, size(hist_dem)[1],T,G)    ######might need to change this later as model changes
    for g in 1:G
        data = hist_dem[:,:,g]
        ## Classify demand and determine forecasting method
        # Outputs idx of demand streams that need to be forecasted by whatever method
        idx_SBA = idClass(data)[7]
        idx_Croston = idClass(data)[8]
        @show(idx_Croston)
        @show(idx_SBA)
        
        ## Forecast SBA and Croston
        # define function handle used to determine alpha value that minimizes mse
        mseh(x,y) = mean(filter(!isnan, (x .- y).^2)) # mse fh. x: historical data, y: forecasted value output.
        # set starting alpha value
        a = [0.1]
        F =  Array{Float64}(undef, (length(idx_SBA) + length(idx_Croston)))
        aoptSBA =Array{Float64}(undef, length(idx_SBA))
        aoptCrost = Array{Float64}(undef, length(idx_Croston))
        
        
        if !isempty(idx_SBA) #checker to ensure idx_SBA has values in it
            for i = idx_SBA # idx of customer demand for product G that requires SBA forecast
                if i != 0 # only do this for nonzero values
                    # Objective is to minimize the mse by optimizing alpha value
                    temp_aoptSBA = Optim.optimize(alpha -> mseh(data[i,:], CrostonSBA(data[i,:], T, alpha)[2]), a, Newton()).minimizer # minimize mse by changing alpha value
                    aoptSBA[i] = temp_aoptSBA[1] # save optimal alpha value for checking
                    F[i,1] = CrostonSBA(data[i,:],T,aoptSBA[i])[1] # Run function using the optimized alpha value
                end
            i # checking
            end
        end
        
        if !isempty(idx_Croston) #checker to ensure idx_Croston has values in it
            for i = idx_Croston # idx of customer demand for product G that requires SBA forecast
                if i != 0 
                    #min mse by changing alpha value
                    temp_aoptCrost = Optim.optimize(alpha -> mseh(data[i,:], Croston(data[i,:],T,alpha)[2]), a, Newton()).minimizer
                    aoptCrost[i] = temp_aoptCrost[1] # save optimal alpha value for checking
                    F[i,1] = Croston(data[i,:],T,aoptCrost[i])[1]
                end
                #i # checking
            end
        end
        #Forecast for time horizon
        forecast[:,:,g] .= ceil.(F)
    end
    
    #now we have the forecast for the DC, we need to find it for each customer
    #This is based on which DC is connected to what customers
    forecast_by_cust = zeros(ncust,T,G)
    
    #building indexes varialbes that is a tuple of length nsup that will hold the indices of which customers each supplier has a valid connection with 
    #this is used to later for a random draw in dishing out demand to customers. 
    indexes = []
    for ns = 1:nsup
        temp = []
        for nc = 1:ncust
            if Suprte[ns,nc] == 1
                push!(temp, nc)
            end
        end
        push!(indexes, temp)
    end
    
    #time to dish out demand to customers
    for t = 1:size(forecast,2) #initiate loop to get through each time period to split up demand appropriately. 
        for g = 1:G #loop thorugh each product
            for ns=1:nsup #loop through each supplier
                if forecast[ns,t,g] > 0
                    n = forecast[ns,t,g] #qty demand to dish out
                    for it = 1:n
                        rand_cust = rand(indexes[ns])
                        forecast_by_cust[rand_cust,t,g] = 1 + forecast_by_cust[rand_cust,t,g]
                    end
                end
            end
        end
    end
    return(forecast_by_cust)
end


"""   **************************************************** THIS FUNCTION IS NOT USED, OLD WAY AND NO WAY TO AGGREGATE DEMAND FOR DCs. LOOK AT NEXT ONE
###### Determining the min policy code
#Markov Bootstrap code for little s inventory policy 
NOTE: this is thought to be from the customers perspective so demand is aggregated for the DCs and the min policy is divided up evenly. This is not used in practice in the current version of the code.  
Inputs:
    nsup: number of suppliers
    ncust: number of customers
    LT: based on LT of resupply from TM to the supplier
    rep: Number of replications to generate LTD based on
    hist_dem: historical demand calculated earlier
Outputs:
    LTD - Lead time duration for each product
    Q - Markov bootstrap for each product 75th percentile  ###### better description maybe needed
    pAll - little s for each product (nsup x G matrix)
"""
function little_s_inventory(nsup, ncust, LT, rep, hist_dem)
    historical_byproduct = sum(hist_dem,dims=1) # forecast based on all customers combined. Look at each product seperateley
    
    #determine LTD (Q) for each product
    G = size(hist_dem)[3]
    
    #initiate arrays
    Q = zeros(G,1)
    LTD = zeros(G,rep)
    for g = 1:G
        (LTD[g,:],Q[g,1]) = MarkovBootstrap(historical_byproduct[:,:,g],LT,rep) #ltd PER PRODUCT
    end
    Q = round.(Q) #this is by product for the entire system. Divide it up between all sups
    
    # setting little s stored as pAll variable
    z = round.(Q/nsup)
    pAll = zeros(nsup,G)
    for g in 1:G
        pAll[:,g] = repeat([z[g]],nsup,1)
    end
    return (LTD, Q, pAll)
end



"""
Determines min policy, MarkovBootstrapping little s. This is from DCs perspective. 
Inputs:
    nsup: number of suppliers
    ncust: number of customers
    LT: based on LT of resupply from TM to the supplier
    rep: Number of replications to generate LTD based on
    hist_dem: historical demand calculated earlier
Outputs:
    LTD - Lead time duration for each product
    Q - Markov bootstrap for each product 75th percentile  ###### better description maybe needed
    pAll - little s for each product (nsup x G matrix)
"""
function little_s_inventory_DC(nsup,ncust,LT,rep,hist_dem,LTreactive=true)
    pAll = zeros(nsup,G)
    Q = zeros(G,1)
    LTD = zeros(G,N,nsup)
    if LTreactive == true
        for ns=1:nsup
            for g = 1:G
                (LTD[g,:,ns],Q[g,1]) = MarkovBootstrap(hist_dem[ns,:,g],LT[1,ns],N)
                Q = round.(Q)
                pAll[ns,g] = Q[g,1]
            end
        end
    else 
        for ns=1:nsup
            for g = 1:G
                (LTD[g,:,ns],Q[g,1]) = MarkovBootstrap(hist_dem[ns,:,g],LT[1],N) #just a sneaky way to make it non-reactive to lead time. Needs to be modified if the first one is meant to be different in an expanded network. 
                Q = round.(Q)
                pAll[ns,g] = Q[g,1]
            end
        end
    end
    return (LTD,Q,pAll)
end


"""
Function that based on number of suppliers and number of customers provided builds Network design where it has 1 for AM capabilities or 0 o/w. Gives all the possible network designs.
"""
function AM_design(nsup,ncust)
    q = []
    for i in 0:2^(nsup+ncust)-1
        push!(q,string(i,base=2,pad=nsup+ncust))
    end
    NetDesign = zeros(2^(nsup+ncust),nsup+ncust)
    for i in 1:length(q)
        for j in 1:length(q[i])
            NetDesign[i,j] = parse(Int8,q[i][j])
        end
    end
    return NetDesign # All suppliers are first, customers are last.
end


"""
Function to split demand randomly among the SLs after assigning the intermittent demand to the DCs.
"""
function split_demand(nreal,T,ncust,nsup,G,Suprte,demand)
    #initiate customer demand array
    cust_dem = zeros(nreal,T,ncust,G)

    #building indexes varialbes that is a tuple of length nsup that will hold the indices of which customers each supplier has a valid connection with 
    #this is used to later for a random draw in dishing out demand to customers. 
    indexes = []
    for ns = 1:nsup
        temp = []
        for nc = 1:ncust
            if Suprte[ns,nc] == 1
                push!(temp, nc)
            end
        end
        push!(indexes, temp)
    end
    
    #for each demand observation here, we randomly dish it out to customers. 
    for real = 1:nreal
        for ns = 1:nsup
            for g=1:G
                for t = 1:size(demand,2)
                    if demand[real,t,ns,g] > 0
                        n = demand[real,t,ns,g] #Demand to dish out to customers
                        for it = 1:n
                            rand_cust = rand(indexes[ns])
                            cust_dem[real,t,rand_cust,g] = 1 + cust_dem[real,t,rand_cust,g]
                        end
                    end
                end
            end
        end
    end
    return cust_dem
end