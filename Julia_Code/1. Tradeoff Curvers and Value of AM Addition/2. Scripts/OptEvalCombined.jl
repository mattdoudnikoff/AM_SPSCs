#Packages needed
using Optim
using Plots
using StatsBase
using Distributions
using Statistics
using Random
using JuMP, Gurobi, GLPK
using TickTock
using Distributed
using CSV
using DataFrames
using JLD2
using Dates
####### This function just checks the OS and helps automate things instead of needed to copy and paste paths
function checkOS()
    path = pwd()
    if path[1] == '/'
        OS = "Mac"
    else
        OS = "Windows"
    end
    return OS
end
OS = checkOS()
file_path = pwd();

#####Load correct .jl files for functions (accommodates for Mac or Windows)
if OS == "Mac"
    inipath = file_path[1:findlast("Julia_Code",file_path)[end]]*"/"
    include(inipath * "1. Tradeoff Curvers and Value of AM Addition/2. Scripts/Generate_Forecast.jl")
    include(inipath * "3. Functions/MILPSetups.jl")
    include(inipath * "3. Functions/Functions.jl")
else 
    inipath = file_path[1:findlast("Julia_Code",file_path)[end]] * "\\"
    include(inipath * "1. Tradeoff Curvers and Value of AM Addition\\2. Scripts\\Generate_Forecast.jl")
    include(inipath * "3. Functions\\MILPSetups.jl")
    include(inipath * "3. Functions\\Functions.jl")
end

function full_run_varied_dem_inputs(nsup, ncust, nTM, valid_TMrte, Suprte, CtTM, CtS, avg_dem, adi, cv2, lev, T, G, obs, nreal, LT, N, AM_capacity, description, NetDes = 0, reactive=true)



    
    #Cost and Lag Inputs (Edit if you really want to)
    ##### NOTE: you can input the exact dimensions and change cost as you please or you can have just a single number there and that will address all the variation in check dims command
    LagTM = LT # Rows = TM Columns = Sup 

    LagS = [2] # Rows = Sup  Columns = Cust 


    CpTM = [1] #rows = num products

    #********CpAM incorporates both supplier and customers becasue AM is located at both locs
    CpAM = [1] #Rows = Sup and Cust, Sups first, Cols = Time, 3D=Prod  

    ##### This might be needed later 
    CtSexpanded = repeat(CtS,1,T,G)

    #All needs to be a function of value of the product. Wont be the same for all products usually
    Ci = [3] #Rows = Sup, Cols = Time, 3D = prod

    CiRetailer = [5] #rows = Retailer, Cols = Time, 3D = prod
    CBO = [378] #Rows = Retailer, Cols = Time, 3D = prod


    ##### Setting up expedited infrastructure  ****** NEEDS WORK TO IMPLEMENT
    CpExp = ones(1,T,G) #covers all expedited prod...... Rows = nTM, Cols = Time, 3D = prod

    ### cost of expedited teransportation. TM direct to SL  ****** NEEDS WORK TO IMPLEMENT
    CtExp = [100] #Transport from RETAILER DIRECT location. Rows = RETAILER Locs, Cols =Time (in the end), 3D Prod (new TM)

    LagExp = [5] #Rows = TM, Cols = Retailers  ****** NEEDS WORK TO IMPLEMENT


    ###### Transhipment  ********* NEEDS WORK TO IMPLEMENT
    #Expedited transport cost from Supplier to Supplier. ONLY USED IN EVAL MILP
    CtransferExp = [0 25 25; 25 0 25; 25 25 0] #square matrix of ALL suppliers. cost of lateral transfer between suppliers

    #valid_transferExprte = [0 1 1; 1 0 1; 1 1 0]
    #valid_transferExprte = [0 0 0; 0 0 0; 0 0 0]

    LagTransferExp = 1 #time lag to transfer between suppliers (ASSUME ALL LAGS ARE 1) ****** NEEDS WORK TO IMPLEMENT

    #The transferred prodcut still flows to the customer in the same manner as normal product. Break it out in separate vars for visualization  ****** NEEDS WORK TO IMPLEMENT
    CtransferExpFLOW = CtSexpanded

    #AM stuff must incorporate both suppliers and customers because both locs have AM, 
    AM_fixed_cost = [121] #suppliers must come first #this is per day

    #Suppliers and cust. Suppliers must come first
    K_stable = [1] #AM indicator. Rows = Sup and Cust. Cols = Time. 3D=Products  ####used to be capacity but i made that a scalar above

    #Expedited Capacity (can change to help prevent eval model from expediting everything)  ****** NEEDS WORK TO IMPLEMENT
    CtExpCapacity = zeros(size(CtExp))
    CtExpCapacity .= 0 #change as needed

    
    
    
    if length(adi) > 1 || length(adi) != length(cv2) || length(cv2) != length(description)
        #generate the initial demand 
        hist_dem = gen_historical(adi[1], cv2[1], lev[1], G, nsup)
        DC_dem, baseline_dem = gen_dem_realization(nreal, adi[1], cv2[1], lev[1], T, G, nsup);
        dem_real = split_demand(nreal,T,ncust,nsup,G,Suprte,DC_dem)
        forecast = forecast_dem(hist_dem,G,T,ncust,Suprte)
        (LTD, Q, pAll) = little_s_inventory_DC(nsup, ncust, LagTM, N, hist_dem,reactive)
        #added to test specific network designs if you know it and want to reduce time/number of runs
        NetDesign = AM_design(nsup,ncust)
        if NetDes == 0
        else
            NetDesign = NetDesign[NetDes,:]
        end
        ndesign = size(NetDesign)[1]
        
        yinit = pAll .+ 10 #Start with more products?
    #LagTM,LagS,CpTM,valid_TMrte,valid_Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable = check_dims(nsup,ncust,nTM,T,G,LagTM,LagS,CpTM,valid_TMrte,Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable)
        
        #Run Opt policy and eval code 
        Big_S,TC_stored = big_S_MILP(T,G,nsup,ncust,CpTM,CpAM,CtTM,CtS,Ci,CiRetailer,CBO,valid_TMrte,Suprte,CtExp,CtransferExp,K_stable,LagS, LagTM,ndesign,AM_fixed_cost,LagTransferExp,CtransferExpFLOW, yinit, pAll, forecast, NetDesign, AM_capacity)

        stored_TC, stored_delivered, stored_BOtotal, stored_BOtrue, stored_BOlag, stored_AMprod, stored_transDC, stored_TMprodCost, 
        stored_TMtransCost, stored_SupTransCost, stored_InvCarCost, stored_AMprodcost, stored_reorder_qty, stored_transDC, stored_supinv, 
        stored_invDC, stored_reorder_delivered  = MILP2(T,G,nsup, ncust,CpTM,CpAM,CtTM,CtS,Ci,CiRetailer,CBO,valid_TMrte,Suprte,CtExp,
        CtransferExp,K_stable,LagS,LagTM,ndesign, AM_fixed_cost,LagTransferExp,CtransferExpFLOW, yinit, pAll, dem_real, NetDesign, 
        AM_capacity, Big_S) 

        TC,BOtot, AMproduced, TMprodCost, TMtransCost, SupTransCost, InvCarCost, AMprodCost, percent_byAM, TCstdev, BOstdev, BO_by_SL, TMcoststdev, TMtransstdev, Suptransstdev, InvCoststdev, AMprodcoststdev = analyze_outputs(stored_TC,stored_delivered,
        stored_BOtrue,stored_AMprod,stored_transDC, stored_TMprodCost, stored_TMtransCost, stored_SupTransCost, stored_InvCarCost, stored_AMprodcost,dem_real)
    
        save_data(description[1],adi[1],cv2[1],avg_dem, TC, BOtot, AMproduced, TMprodCost, TMtransCost, SupTransCost, InvCarCost, AMprodCost, Big_S, pAll, percent_byAM, TCstdev, BOstdev, BO_by_SL, TMcoststdev, TMtransstdev, Suptransstdev, InvCoststdev, AMprodcoststdev)
        
        for i = 2:length(adi)
            #generate the  demand 
            hist_dem = gen_historical(adi[i], cv2[i], lev[i], G, nsup)
            DC_dem, extra = gen_dem_realization(nreal,adi[i],cv2[i],lev[i],T,G,nsup,baseline_dem) #THIS IS THE CHANGE TO KEEP DEMAND LEVELS WITHIN 5%
            dem_real = split_demand(nreal,T,ncust,nsup,G,Suprte,DC_dem)
            forecast = forecast_dem(hist_dem,G,T,ncust,Suprte)
            (LTD, Q, pAll) = little_s_inventory_DC(nsup, ncust, LagTM, N, hist_dem,reactive)
            
            yinit = pAll #Start with more products?
    #LagTM,LagS,CpTM,valid_TMrte,valid_Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable = check_dims(nsup,ncust,nTM,T,G,LagTM,LagS,CpTM,valid_TMrte,Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable)

            #Run Opt policy and eval code 
            Big_S,TC_stored = big_S_MILP(T,G,nsup,ncust,CpTM,CpAM,CtTM,CtS,Ci,CiRetailer,CBO,valid_TMrte,Suprte,CtExp,CtransferExp,K_stable,LagS, LagTM,ndesign,AM_fixed_cost,LagTransferExp,CtransferExpFLOW, yinit, pAll, forecast, NetDesign, AM_capacity)

            stored_TC, stored_delivered, stored_BOtotal, stored_BOtrue, stored_BOlag, stored_AMprod, stored_transDC, stored_TMprodCost, 
            stored_TMtransCost, stored_SupTransCost, stored_InvCarCost, stored_AMprodcost, stored_reorder_qty, stored_transDC, stored_supinv, 
            stored_invDC, stored_reorder_delivered  = MILP2(T,G,nsup, ncust,CpTM,CpAM,CtTM,CtS,Ci,CiRetailer,CBO,valid_TMrte,Suprte,CtExp,
            CtransferExp,K_stable,LagS,LagTM,ndesign, AM_fixed_cost,LagTransferExp,CtransferExpFLOW, yinit, pAll, dem_real, NetDesign, 
            AM_capacity, Big_S) 

            TC,BOtot, AMproduced, TMprodCost, TMtransCost, SupTransCost, InvCarCost, AMprodCost, percent_byAM, TCstdev, BOstdev, BO_by_SL, TMcoststdev, TMtransstdev, Suptransstdev, InvCoststdev, AMprodcoststdev = analyze_outputs(stored_TC,stored_delivered,
            stored_BOtrue,stored_AMprod,stored_transDC, stored_TMprodCost, stored_TMtransCost, stored_SupTransCost, stored_InvCarCost, stored_AMprodcost,dem_real)

            save_data(description[i],adi[i],cv2[i],avg_dem, TC, BOtot, AMproduced, TMprodCost, TMtransCost, SupTransCost, InvCarCost, AMprodCost, Big_S, pAll, percent_byAM, TCstdev, BOstdev, BO_by_SL, TMcoststdev, TMtransstdev, Suptransstdev, InvCoststdev, AMprodcoststdev)
        end
    else
        error("ERROR: Length of ADI is not the same as CV2 or description. Edit Accordingly")
    end
end