"""
    This function takes inputs, checks the dimensions and make sure they are right for the MILP.
"""
function check_dims(nsup,ncust,nTM,T,G,LagTM,LagS,CpTM,valid_TMrte,Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable)
    if size(LagTM) != (nTM, nsup)
        if length(LagTM) == 1
            LagTM = repeat(LagTM, nTM, nsup)
        else
            error("ERROR: LagTM length is not 1 or the size nTM x nsup. Edit Accordingly")
        end
    end
    
    if size(LagS) != (nsup,ncust)
        if length(LagS) == 1
            LagS = repeat(LagS, nsup,ncust)
        else
            error("ERROR: LagS length is not 1 or the size of nsup x ncust. Edit Accordingly")
        end
    end
    
    if length(CpTM) == 1
        CpTM = repeat(CpTM, G,T)    ######might have to look at this later
    elseif length(CpTM) == nsup
        CpTM = repeat(CpTM',1,T) ####don't think this actually needs to change because we assume constant
    else
        error("ERROR: CpTM length is not 1 or the size of nsup. Edit Accordingly")
    end
    
    
    #### Maybe change this later but need to change the MILP too 
    valid_TMrte = repeat(valid_TMrte .* Inf, 1, T, 1)
    valid_TMrte[isnan.(valid_TMrte)] .= 0
    
    valid_Suprte = zeros(nsup,T,G,ncust) #can drop the T dimension and make MILP simpler
    
    for i = 1:nsup
        for j=1:ncust
            valid_Suprte[i,:,:,j] .= (Suprte[i,j] * Inf)
        end
    end
    valid_Suprte[isnan.(valid_Suprte)] .= 0 #replace NaNs with 0
    
    if length(Ci) == 1
        Ci = repeat(Ci,nsup,T+1,G)
        Ci[:,1,:] .= 0
    elseif length(Ci) == nsup
        Ci = repeat(Ci,1,T+1,G)
        Ci[:,1,:] .= 0
    else
        error("ERROR: Ci length is not 1 or the size of nsup. Edit Accordingly")
    end
    
    if length(CiRetailer) == 1
        CiRetailer = repeat(CiRetailer,ncust,T+1,G)
        CiRetailer[:,1,:] .= 0
    elseif length(CiRetailer) == ncust
        CiRetailer = repeat(CiRetailer,1,T+1,G)
        CiRetailer[:,1,:] .= 0
    else
        error("ERROR: CiRetailer length is not 1 or the size of ncust. Edit Accordingly")
    end

    if length(CBO) == 1
        CBO = repeat(CBO,ncust,T+1,G) #might be able to drop the T dimensions
    elseif length(CBO) == ncust
        CBO = repeat(CBO,1,T+1,G) #might be able to drop the T dimensions
    else
        error("ERROR: CBO length is not 1 or the size of ncust. Edit Accordingly")
    end
    
    if length(CpAM) == 1
        CpAM = repeat(CpAM,(nsup+ncust),T,G) #might be able to drop the T dimension
    elseif length(CpAM) == (nsup+ncust)
        CpAM = repeat(CpAM,1,T,G) #might be able to drop the T dimension
    else
        error("ERROR: CpAM length is not 1 or the length of nsup+ncust. Edit Accordingly")
    end
    
    if length(AM_fixed_cost) == 1
        AM_fixed_cost = repeat(AM_fixed_cost,ncust+nsup,1)
    elseif length(AM_fixed_cost) == (ncust + nsup)
    else
        error("ERROR: AM_fixed_cost length is not 1 or the length of nsup+ncust. Edit Accordingly")
    end
    
    if length(K_stable) == 1 
        K_stable = repeat(K_stable,ncust+nsup)
    elseif length(K_stable) == ncust+nsup
    else
        error("ERROR: K_stable length is not 1 or length of nsup+ncust. Edit Accordingly.")
    end
    
    return LagTM,LagS,CpTM,valid_TMrte,valid_Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable
end





"""   First MILP
This function runs the first MILP needed to solve for Big S Policy. 

Coefficient Description:
Inputs:
    T
    G
    nsup
    ncust
    CpTM
    CpAM
    CtTM
    CtS
    Ci
    CiRetailer
    CBO
    valid_TMrte
    Suprte
    CtExp
    CtransferExp
    K_stable
    LagS
    LagTM
    ndesign
    AM_fixed_Cost
    LagTransferExp
    CTransferExpFLOW
    yinit
    pAll
    forecast
    NetDesign
    AM_capacity

Used in Obj Function
    CpTM = Cost of Production Tradition Manufacturing (TM) - NumTM x T x G matrix
    CpAM = Cost of Production Additive Manufacturing (AM) - NumPossAMlocs x T x G matrix
    CtTM = Cost of Transportation from TM to Supplier (DCs) - NumTM2DCroutes x T x G matrix 
    CtS = Cost of Transportation from Supplier to SLs - NumDCs x (ncust*T) x G matrix
    Ci = Inventory Cost (at supplier) - nsup x T+1 x G matrix
    CiRetailer = Inventory Cost (at SLs) - ncsut x T+1 x G matrix
    CBO = Backorder Costs - ncust x T+1 x G matrix
Other Coefs
    Creorder - Cost of reordering (not sure if this is actually used)
    CpExp, CtExp = Expedited production + transportation costs used in Eval (TM to SL)
"""
function big_S_MILP(T,G,nsup,ncust,CpTM,CpAM,CtTM,CtS,Ci,CiRetailer,CBO,valid_TMrte,Suprte,CtExp,CtransferExp,K_stable,LagS,LagTM,ndesign,AM_fixed_cost,LagTransferExp,CtransferExpFLOW, yinit, pAll, forecast, NetDesign, AM_capacity)
    
    LagTM,LagS,CpTM,valid_TMrte,valid_Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable = check_dims(nsup,ncust,nTM,T,G,LagTM,LagS,CpTM,valid_TMrte,Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable)
    
    #CtExp = repeat(CtExp',1,T,G)  #Columns turn to time
    
    #CtransferExp = repeat(CtransferExp, 1, T, G) #Def needs to be modified later. Haven't tested these kinks out.
    
    ###### MIGHT BE NEEDED FOR TRANSHIPMENT BUT NEED TO WORK ON
    #valid_transferExprte = repeat(valid_transferExprte*Inf,1,T,G)
    #valid_transferExprte[isnan.(valid_transferExprte)].=0 #replace NaNs with 0
    
    #Lag Multiplier (to account for supplier lag in shipping to SLs, not a problem for DCs as current problem stands, because of one TM)
    LagSMult = ones(nsup,T,G,ncust)
    for i=1:nsup
        for j=1:ncust
            for g=1:G
                for t=1:LagS[i,j]
                    LagSMult[i,t,g,j] = 0
                end
            end
        end
    end
    K_stable = repeat(K_stable,1,T,G)

    MinPolicy = pAll #Min policy (little s)

    lbinvenDC = zeros(nsup,T+1,G)
    ubinvenDC = zeros(nsup,T+1,G) .+ (1*Inf)
    lbinvenSL = zeros(ncust,T+1,G)
    ubinvenSL = zeros(ncust, T+1,G) #.+ (Inf)
    for g=1:G
        for i=1:nsup
            lbinvenDC[i,1,g] = yinit[i,g] 
            lbinvenDC[i,T+1,g] = yinit[i,g]
            ubinvenDC[i,1,g] = yinit[i,g]
            ubinvenDC[i,T+1,g] = yinit[i,g]
        end
        #### Had this loop before but it breaks it when I changed to DC perspective.
        #for i=1:ncust
            #lbinvenSL[i,1,g] = forecast[i,1,g]
            #ubinvenSL[i,1,g] = forecast[i,1,g]
        #end
    end
    
    TC_stored = zeros(ndesign)
    Big_S = zeros(nsup,G,ndesign)
    # Iterate through network designs and save the policies at the end
    for d=1:ndesign
        println("Design ", string(d))
        #println("Check Start Loop")
        AM_force = NetDesign[d,:] #Changes a network design each time (AM locations)
        K = K_stable .* AM_force # Eliminates capacity of non-AM locs. (Turn AM on/off)
        #Keeps the AM capacity of the system constrant (Excluding the No AM case)
        if d == 1
            AM_cap = AM_force
        else
            AM_cap = AM_force .* floor.(AM_capacity ./ sum(AM_force))
        end
        
        model1 = Model(Gurobi.Optimizer)
        #model1 = Model(GLPK.Optimizer)
        #Setting up sets for calcualtions
        NS = 1:nsup
        NC = 1:ncust
        I = 1:(nsup+ncust)
        ## Set up Decision Variables
        #println("Check DVs")
        @variable(model1, 0 <= TMprod[g=1:G,t=1:T],Int) #quantity of part g traditionally produced in period t at loc i in P
        @variable(model1, 0 <= AMprod[i=1:(nsup+ncust),t=1:T,g=1:G],Int) #quantity of part g produced with AM in period t at location i at DC or Service Loc
        @variable(model1, 0 <= transTM[ns=1:nsup,t=1:T,g=1:G] <= valid_TMrte[ns,t,g],Int) #flow of part g transported from location TM to DC in time t
        @variable(model1, 0 <= transDC[ns=1:nsup,t=1:T,g=1:G,nc=1:ncust] <= valid_Suprte[ns,t,g,nc],Int) #flow of part g from m DC loc to i SL loc at time t.
        @variable(model1, lbinvenDC[ns,t,g] <= invenDC[ns=1:nsup,t=1:T+1,g=1:G] <= ubinvenDC[ns,t,g],Int) #quantity of part g carried as inventory at location i in DC from period t-1 to t
        @variable(model1, lbinvenSL[nc,t,g] <= invenSL[nc=1:ncust,t=1:T+1,g=1:G] <= ubinvenSL[nc,t,g],Int) #quantity of part g carried as inventory at location i at SLs from period t
        @variable(model1, 0 <= backorder[nc=1:ncust,t=1:T+1,g=1:G],Int) #quantity of unsatisfied demand of part g at location i in SL from period t
        @variable(model1, 0 <= reorderQTY[ns=1:nsup,t=1:T,g=1:G],Int) #quantity to reorder for supplier i at time t for product k
        @variable(model1, 0 <= delivered[nc=1:ncust,t=1:T,g=1:G],Int) #quantity of item g "delivered/used" at SL i at time t
        @variable(model1, p[ns=1:nsup,t=1:T,g=1:G], Bin) # reordrer indicator for part g at location i (DC) at time period t
        @variable(model1, e[i=1:(ncust+nsup),t=1:T,g=1:G], Bin) #AM production indicator for part g at location i (DCs or SLs) in time period t
        @variable(model1, MinPolicy[ns,g] <= MaxPolicy[ns=1:nsup, g=1:G] <= Inf, Int) ###############Solving for this, Big S
        #Set up obejctive function
        @objective(model1, Min, 
            sum(CpTM[g,t]*TMprod[g,t] for g=1:G,t=1:T) + #TM production cost
            sum(CtTM[ns]*transTM[ns,t,g] for ns=1:nsup, t=1:T, g=1:G) +  #transportation cost TM to supplier
            sum(CtS[ns,nc]*transDC[ns,t,g,nc] for ns=1:nsup, t=1:T, g=1:G, nc=1:ncust) + #transportation supplier to SLs
            sum(Ci[ns,t,g]*invenDC[ns,t,g] for ns=1:nsup, t=2:T+1, g=1:G) + #DCs inventory cost
            sum(CiRetailer[nc,t,g]*invenSL[nc,t,g] for nc=1:ncust, t=2:T+1, g=1:G) + #Customer inventory cost
            sum(CpAM[i,t,g]*AMprod[i,t,g] for i=1:(nsup+ncust), t=1:T, g=1:G) + #AM production cost
            sum(CBO[nc,t,g]*backorder[nc,t,g] for nc=1:ncust, t=1:T+1, g=1:G)); #backorder cost
        
        #println("Start Constraints")
        ###Constraints
        #Flow balance produced at time T = ∑ of all things transported at time T
        #i.e. No inventory stored at TM
        for t = 1:T
            for g = 1:G
                @constraint(model1, TMprod[g,t] == sum(transTM[ns,t,g] for ns ∈ NS)) 
            end
        end

        #Inventory Balance at DCs
        for i=1:nsup
            for g=1:G
                for t=1:LagTM[i]
                    @constraint(model1, invenDC[i,t+1,g] == invenDC[i,t,g] + AMprod[i,t,g] - sum(transDC[i,t,g,nc] for nc ∈ NC))
                end
                for t=(LagTM[i]+1):T
                    @constraint(model1, invenDC[i,t+1,g] == invenDC[i,t,g] + AMprod[i,t,g] + sum(transTM[i,t-LagTM[i],g]) - sum(transDC[i,t,g,nc] for nc ∈ NC))
                end
            end
        end
        
        
        #Inventory Balance at SLs
        for i=1:ncust
            for g=1:G
                for t=1:T
                    @constraint(model1,invenSL[i,t+1,g] == invenSL[i,t,g] + AMprod[i+nsup,t,g] + 
                        sum(LagSMult[ns,t,g,i] * transDC[ns,maximum([t-LagS[ns,i],1]),g,i] for ns ∈ NS) - delivered[i,t,g])
                end
            end
        end

        #Backorder
        for i=1:ncust
            for t=1:T
                for g=1:G
                    @constraint(model1,backorder[i,t+1,g] == backorder[i,t,g] + forecast[i,t,g] - delivered[i,t,g])
                end
            end
        end

        #Trigger re-order constraints
        for i=1:nsup
            for g=1:G
                @constraint(model1, MinPolicy[i,g] <= invenDC[i,2,g] + (100000 * p[i,1,g]))
                @constraint(model1, MinPolicy[i,g] >= invenDC[i,2,g] - (100000 * (1-p[i,1,g])))
                for t=2:T
                    @constraint(model1, MinPolicy[i,g] <= invenDC[i,t+1,g] + (100000 * (p[i,t,g])) + sum(reorderQTY[i,m,g] for m=(t-1):maximum([t-LagTM[i],1])))
                    @constraint(model1, MinPolicy[i,g] >= invenDC[i,t+1,g] - (100000 * (1-p[i,t,g])) + sum(reorderQTY[i,m,g] for m=(t-1):maximum([t-LagTM[i],1])))
                end
            end
        end

        ##### Reorder Constraints
        for i=1:nsup
            for t=1:T
                for g=1:G
                    #Reorder quantity
                    @constraint(model1, reorderQTY[i,t,g] == (MaxPolicy[i,g] - invenDC[i,t+1,g]) * p[i,t,g])

                    #Replenishment constraints
                    @constraint(model1, reorderQTY[i,t,g] <= 100000 * p[i,t,g])
                    @constraint(model1, reorderQTY[i,t,g] <= MaxPolicy[i,g] - invenDC[i,t+1,g])
                    @constraint(model1, reorderQTY[i,t,g] >= MaxPolicy[i,g] - invenDC[i,t+1,g] - 10000*(1-p[i,t,g]))
                end
            end
        end

        #Linking Reorder to trans from TM to DCs
        for i=1:nsup
            for t=1:T
                for g=1:G
                    @constraint(model1, transTM[i,t,g] == reorderQTY[i,t,g])
                end
            end
        end

        #AM capacity and production indicator
        for i=1:(nsup+ncust)
            for t=1:T
                for g=1:G
                    @constraint(model1, AMprod[i,t,g] <= K[i,t,g] * AM_cap[i] * e[i,t,g]) #This constraint links production indicator to AM production
                end
                #link e with as a production indicator (wheather or not AM is currently being used for a product)
                @constraint(model1, sum(e[i,t,g] for g=1:G) == 1)
            end
        end
        #AM capacity constraint
        for t=1:T
            @constraint(model1, sum(AMprod[i,t,g] for i=1:(nsup+ncust), g=1:G) <= AM_capacity)
        end
        set_optimizer_attribute(model1, "OutputFlag", 1)
        set_optimizer_attribute(model1, "MIPGap", .05) #5% optimality gap
        #set_optimizer_attribute(model1, "TimeLimit", )
        set_silent(model1) #can take this out, but runs faster when you don't print the messages/results
        optimize!(model1);
        #println("Check Done")
        TC_stored[d] = objective_value(model1)

        Big_S[:,:,d] = value.(MaxPolicy)
    end
    
    return Big_S, TC_stored
end
    













""" Second MILP
    This function runs the second MILP that looks at day by day functions and makes a decision everyday then eventually calculates many Costs and returns them. 
"""
function MILP2(T,G,nsup,ncust,CpTM,CpAM,CtTM,CtS,Ci,CiRetailer,CBO,valid_TMrte,Suprte,CtExp,CtransferExp,K_stable,LagS,LagTM,ndesign,AM_fixed_cost,LagTransferExp,CtransferExpFLOW, yinit, pAll, dem_real, NetDesign, AM_capacity, Big_S) 
    #### Eval Code
    #nd = 16  #Use this if testing single runs to check for speed 
    #initialize all stored variables for later
    maxT = maximum([maximum(LagS) maximum(LagTM)]) #some arrays need to be past time horizon because of shipping lag, this is the maxT to save space later
    stored_transDC = zeros(nsup,T,G,ncust,nreal,ndesign);
    stored_AMprod = zeros(nsup+ncust,T,G,nreal,ndesign);
    stored_SupTransCost = zeros(T,nreal,ndesign);
    stored_AMprodcost = zeros(T,nreal,ndesign);
    stored_TMprodCost = zeros(T,nreal,ndesign);
    stored_TMprod = zeros(T,nreal,ndesign);
    stored_InvCarCost = zeros(T,nreal,ndesign);
    stored_TC = zeros(T,nreal,ndesign);
    stored_reorder_qty = zeros(nsup,T,G,nreal,ndesign);
    stored_TMtransCost_Itemized = zeros(nsup,T,nreal,ndesign);
    stored_delivered = zeros(ncust,T+maxT,G,nreal,ndesign);  #####may need to change the max later
    #stored_reorder = zeros();
    stored_reorder_delivered = zeros(nsup,T+maxT,G,nreal,ndesign); 
    stored_TMtransCost = zeros(T,nreal,ndesign);
    stored_BOtotal = zeros(ncust,T+1,G,nreal,ndesign);
    stored_BOtrue = zeros(ncust, T+1, G, nreal,ndesign);
    stored_BOlag = zeros(ncust,T+1,G,nreal,ndesign);
    stored_supinv = zeros(nsup,T+maxT,G,nreal,ndesign);
    stored_invDC = zeros(nsup,T+maxT,G,nreal,ndesign);


    LagTM,LagS,CpTM,valid_TMrte,valid_Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable = check_dims(nsup,ncust,nTM,T,G,LagTM,LagS,CpTM,valid_TMrte,Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable)
    
    #### if not solving take "@distributed" out and error should pop up. 
    for nd = 1:ndesign
        #update AM force
        println("Design ", string(nd))
        AM_force = NetDesign[nd,:] #changes the network design of AM lcoations
        
        #Keeps the AM capacity of the system constrant (Excluding the No AM case)
        if nd == 1
            AM_cap = AM_force
        else
            AM_cap = AM_force .* floor.(AM_capacity ./ sum(AM_force))
        end

        #Setting Min/Max Policy found in MILP
        MaxPolicy = Big_S[:,:,nd]
        if nsup == 1
            MaxPolicy_3d = reshape(MaxPolicy,(nsup,1,G)) #changes to be rows=nsup, cols=t, 3d=product
            MinPolicy = pAll
            MinPolicy_3d = reshape(MinPolicy,(nsup,1,G))
        else
            MaxPolicy_3d = reshape(MaxPolicy,(nsup,1,G)) #changes to be rows=nsup, cols=t, 3d=product
            MinPolicy = pAll
            MinPolicy_3d = reshape(MinPolicy,(nsup,1,G))
        end
            

        #Loop for demand realizations
        for dr = 1:nreal
            dem_running = dem_real[dr,:,:,:]
            #needs to be reset every run 
            #starting with 5 more products than little s
            yinit = MinPolicy .+ 5
            
            #setting up inventories
            sup = zeros(nsup,T+maxT,G)
            #setup so sup is nsup, time, by nprods matrix,
            for i = 1:nsup
                for g = 1:G
                    sup[i,1,g] = yinit[i,g]
                    stored_supinv[i,1,g,dr,nd] = yinit[i,g]
                end
            end


            previous_orders = zeros(nsup,T+maxT,G) # reset for each realization
            reorder_delivered = zeros(nsup,T+maxT,G) # Destination  matrix for when the re-ordered product will be delivered. Rows = Suppliers. Columns = Time

            #iterate through time
            for t=1:T
                #println("This is time period ", string(t))
                ## Reduce dimensions of select variables for evaluation in MILP
                # Assumption CtS, CBO, CpAM, K do not change over time. If they do then
                # just need to extract the appropriate columns on the appropriate loop
                # through time.
                dem = dem_running[t,:,:];# extract the current time period for use in MILP
                Ci2 = Ci[:,t:t+1,:]; # reduce to 2
                Ci2[:,1,:] .= 0 # First time period (yinit = 0 cost)
                CBO2 = CBO[:,t:t+1,:]; # reduce to 2
                CBO2[:,1,:] .= 0;

                CpAM2 = CpAM[:,t,:] # reduce to 1

                # Reduce to 1
                # Eliminate Capacity of Non-AM 'forced' Sites
                K = K_stable.*AM_force;


                CiRetailer2 = CiRetailer[:,t:t+1,:]; # reduce to 2
                CiRetailer2[:,1,:] .= 0;



                # Not used in MILP but used in future calcs
                CpTM2 = CpTM[:,t,:];

                # Reduce Expedited TM-Cust
                #CpExp2 = CpExp[:,t,:]; # reduce to 1
                #CtExp2 = CtExp[:,t,:]; # reduce to 1

                #CtExpCapacity2 = CtExpCapacity[:,t,:];


                # Reduce the valid rte UBs for use in MILP
                valid_Suprte2 = valid_Suprte[:,t,:,:]; # Used in MILP
                valid_TMrte2 = valid_TMrte[:,t,:]; # Not used in this code..


                # Reduce Sup-Sup Exp transfer variables and UBs
                #CtransferExp2 = CtransferExp[:,((t-1)*3+1):(t-1)*3+3,:];
                #valid_transferExprte2 = valid_transferExprte[:,((t-1)*3+1):(t-1)*3+3,:];

                #CtransferExpFLOW2 = CtransferExpFLOW[:,t:t+1,:];

                yinit = sup[:,t,:] # update inventory for current t
                lbSupply = zeros(nsup,2,G)
                ubSupply = ones(nsup,2,G) .* Inf
                for i = 1:nsup
                    for g=1:G
                        lbSupply[i,1,g] = yinit[i,g]
                        ubSupply[i,1,g] = yinit[i,g]
                    end
                end
                #yinit = permute(yinit,[2 1 3]); 
                y = zeros(nsup,T,G)

                #Set up sets (not necessarily needed but is used)
                NC = 1:ncust
                NS = 1:nsup
                
                #println("CHECK FOR YINIT AND LBSUP AND UBSUP")
                #println(yinit)
                #println(lbSupply)
                #println(ubSupply)
                
                #model = Model(Gurobi.Optimizer);
                model = Model(GLPK.Optimizer)
                ## Set up Decision Variables
                 #mp = buildMILP(Ci,CpAM,CtS,cdel,CiRetailer,CBO,e,CpExp,CtExp,CtransferExp,CtransferExpFLOW)
                @variable(model, lbSupply[ns,tt,g] <= invenDC[ns=1:nsup,tt=1:2,g=1:G] <= ubSupply[ns,tt,g]) #inventory quantity of prod g at DC nc at time t
                @variable(model, 0 <= AMprod[i=1:(nsup+ncust),g=1:G] <= AM_cap[i]) #quantity of product g produced at location i
                #@variable(model, 0 <= transTM[i=1:nsup,g=1:G]) 
                @variable(model, 0 <= transDC[ns=1:nsup,g=1:G,nc=1:ncust] <= valid_Suprte2[ns,g,nc]) #quantity of product g transported from DC ns to SL nc
                @variable(model, 0 <= delivered[nc=1:ncust, g=1:G]) #Quantity of product g delivered/utilized at SL nc
                @variable(model, 0 <= invenSL[nc=1:ncust,tt=1:2,g=1:G] <= 0) #Inventory at SL
                @variable(model, 0 <= backorder[nc=1:ncust,tt=1:2,g=1:G]) #Backorder DV 
                @variable(model, e[i=1:(ncust+nsup),g=1:G], Bin) #AM production indicator
                #@variable(model, 0 <= ExpProd[i=1:ncust,g=1:G]) #Expedited Production
                #@variable(model, 0 <= transExp[i=1:ncust,g=1:G]) #Expedited transportation
                #variable(model, transferExp[]) #Lateral Transfer between DCs

                #Set up obejctive function
                @objective(model, Min, 
                    sum(Ci2[ns,tt,g] * invenDC[ns,tt,g] for ns=1:nsup, tt=1:2, g=1:G) + #DC inventory cost
                    sum(CpAM2[i,g] * AMprod[i,g] for i=1:(nsup+ncust), g=1:G) + #AM production cost
                    sum(CtS[ns,nc] * transDC[ns,g,nc] for ns=1:nsup, g=1:G, nc=1:ncust) + #DC to SL transportation cost
                    sum(CiRetailer2[nc,tt,g] * invenSL[nc,tt,g] for nc=1:ncust, tt=1:2, g=1:G) + #SL inventory cost
                    sum(CBO2[nc,tt,g] * backorder[nc,tt,g] for nc=1:ncust, tt=1:2, g=1:G)) #+ #Backorder Cost
                    #sum(CpExp2[nc,g] * ExpProd[nc,g] for nc=1:ncust, g=1:G) + #Expedited production cost
                    #sum(CtransferExp2[nc,g] * transExp[nc,g] for nc=1:ncust, g=1:G) + #Expedtied transportation cost


                for g = 1:G
                    for tt=1:1
                        # Expedited TM Flow Balance
                        # Amount Exp Produced - sum of Amount EXPEDITED out to all RETAILERS.
                        #@constraint(model, ExpProd[g] == sum(transEXP[i,g] for i=1:nsup))

                        for ns = 1:nsup

                            # Distribution Center (Supplier) Flow Balance
                            # In from TM, + Current Inv + AM Produced - Leftover Inv - Sum of everything shipped out. 
                            @constraint(model, invenDC[ns,tt,g] + AMprod[ns,g] - invenDC[ns,tt+1,g] - sum(transDC[ns,g,nc] for nc ∈ NC) == 0 )
                            # *** want to add: in all in 1st column, all rows and all out: 1st row all columns


                            # Build the transfer connections for 'In and out'. It is OK if the diagnoal is 'wrong' bc the UB on it is zero so nothing will ever flow over it.       
                            #zz = zeros(size(CtransferExp));
                            #zz[:,ns,g] = 1;
                            #zz[ns,:,g] = -1;
                            #mp.addcstr({[1 -1],{ns,[tt tt+1],g}},{ns,tt,g},{-1,{ns,':',g}},0,0,0,0,0,0,{zz,{':',':',g}},{-1,{ns,':',g}},'=',0);

                            #mp.addcstr(0,0,0,0,0,0,0,0,0,{':',ns,g},'=',{ns,':',g}); # Everything in from a lateral transfer must be transported across the special transfer arcs.


                            end #end ns

                        for nc = 1:ncust;
                            # Retailer Flow Balance
                            # Sum of everything in from all suppliers - delivered +Inv In - Inv out
                                # This places a '1' in the appropraite time period that product left the supplier based on the lag.
                            @constraint(model, AMprod[nc+nsup,g] + sum(transDC[ns,g,nc] for ns ∈ NS) - delivered[nc,g] + invenSL[nc,tt,g] - invenSL[nc,tt+1,g] == 0)

                            #BO
                            # in from delivered, -BO in + BO out, = Demand
                            @constraint(model, delivered[nc,g] - backorder[nc,tt,g] + backorder[nc,tt+1,g]== dem[nc, g])
                            end #end nc

                        # Production indicator constraint enforecement
                        #  Fixed Cost no longer active in MILP.
                        # ****************indedxing is different becuase it is for both
                        # suppliers and customers.**********
                        for i = 1:nsup+ncust;
                            @constraint(model, AMprod[i,g] <= K[i] * e[i,g] * AM_cap[i])
                        end
                    end #end of tt loop
                end #end of G loop
                #link e with as a production indicator (wheather or not AM is currently being used for a product)
                for i = 1:(nsup+ncust)
                    @constraint(model, sum(e[i,g] for g=1:G) <= 1)
                end
                #Add total network capacity limit
                @constraint(model, sum(AMprod[i,g] for i=1:(nsup+ncust),g=1:G) <= AM_capacity)
                set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_ALL)
                #set_optimizer_attribute(model, "OutputFlag", 1)
                #set_optimizer_attribute(model, "MIPGap", .001)
                #set_optimizer_attribute(model, "TimeLimit", )
                set_silent(model) #can take this out, but runs faster when you don't print the messages/results
                optimize!(model);

                #store the data 
                stored_transDC[:,t,:,:,dr,nd] .= value.(transDC)
                stored_AMprod[:,t,:,dr,nd] .= value.(AMprod)
                stored_SupTransCost[t,dr,nd] = sum(CtS[ns,nc] * value.(transDC)[ns,g,nc] for ns=1:nsup, g=1:G, nc=1:ncust)
                stored_AMprodcost[t,dr,nd] = sum(CpAM2[i,g] * value.(AMprod)[i,g] for i=1:(nsup+ncust), g=1:G)
                stored_InvCarCost[t,dr,nd] = sum(Ci2[ns,tt,g] * value.(invenDC)[ns,tt,g] for ns=1:nsup, tt=1:2, g=1:G)
                stored_invDC[:,t+1,:,dr,nd] = value.(invenDC)[:,2,:]
                
                #println("Value of invenDC")
                #println(value.(invenDC))
                #println(size(value.(invenDC)))
                

                
                #Update inventory for next time period
                sup[:,t+1,:] = value.(invenDC)[:,2,:] + sup[:,t+1,:] # Represents the amount of inventory on hand to START that period.  Row = Suppliers. Columns = Time. 
                stored_supinv[:,t+1,:,dr,nd] = sup[:,t+1,:]
                ######MIGHT BE OF INTEREST LATER NEED TO BUILD STRUCTURE FOR THIS STORAGE
                ## AM locations used
                #AM_loc_used(:,t,:) = x.CpAM > 0;
                ## AM locations used
                #AM_qty_used(:,t,:) = x.CpAM;

                #Inventory Carried
                
                
                y[:,t,:] = value.(invenDC)[:,2,:] # Inventory carried between periods. (t1 is the inventory that was carried from t1 to t2)
                                # This is different than 'sup' at the end of the loop becuase y extracts your inventory before re-orders arrive.
                                # Use y for inventory carried between periods.
                                # Use sup to see what the actual inventory is at the start of a period (it adds in re-orders)

                # Determine the re-order quantity and place the re-ordered quantiy
                # in sup matrix based on the transportation lag from TM to sup.

                # This should prevent re-order every period, by looking back over
                # the lag time. If your inventory + whatever you have enroute is <
                # Minpolicy then order again
                
                log = zeros(nsup,1,G)
                if t==1
                    log = reshape((y[:,t,:] .+ previous_orders[:,t,:] .< MinPolicy),(nsup,1,G))
                    #println("Inv and ordered")
                    #println((y[:,t,:] .+ previous_orders[:,t,:]))
                else 
                    log = reshape((y[:,t,:] .+ permutedims(sum(previous_orders[:,maximum([1, t-6]):t-1,:],dims=2),[1,3,2])) .< MinPolicy,(nsup,1,G))
                    #println("Inv and order")
                    #println(y[:,t,:] .+ permutedims(sum(previous_orders[:,maximum([1, t-6]):t-1,:],dims=2),[1,3,2]))
                end
                ################### if you change LagTM to non-uniform look at if statement above and edit accordingly (change the 6 somehow), most likely need a loop for each product (if change lag)
 
                
                b = log .* MaxPolicy_3d; # Zeros out corresponding policies that dont need to order
                bb = log .* reshape(y[:,t,:],(nsup,1,G)); # Zeros out corresponding suppliers that dont need to order
                reorder_qty  = b.-bb; # subtract current inventory from policies. This is the re-order qty per supplier.
                previous_orders[:,t,:] = reorder_qty
                
                
                # Now place in the appropriate place (using nested loop to get it right) ##COMMENT: g loop might not be needed because same lag time
                for ns=1:nsup
                    for g=1:G
                        sup[ns,t+LagTM[ns],g] = reorder_qty[ns,1,g] + sup[ns,t+LagTM[ns],g] # Add the re-ordered qty to the location where it will be delivered.
                        stored_supinv[ns,t+LagTM[ns],g,dr,nd] = sup[ns,t+LagTM[ns],g]
                        reorder_delivered[ns,t+LagTM[ns],:] = reorder_qty[ns,1,:] # Just a data tracking matrix. Displays when re-ordered goods, for  each supplier  are  delivered.
                        stored_reorder_delivered[ns,t+LagTM[ns],:,dr,nd] = reorder_qty[ns,1,:]
                    end
                end

                # Save Data  
                stored_reorder_qty[:,t,:,dr,nd] = reorder_qty ; # Reordered quantity per customer per product
                

                """ NOT SURE IF NEEDED

                # NOTE: Sup matrix recieved the re-ordered goods.  At this point,     ################ WORK ON
                # it represents the full amount of inventory available at
                # suppliers for t+1 and on.
                time(t).sup = sup; # This is the inventory available at the beginiing of the next time period time period.
                              # ex. time(1).sup represents the supply conditions to start t2.
                              # NOTE: AM Capacity is added in at the beginning
                              # of next period.
                """

                #Costs associated with re-orders
                
                # Production Cost TM
                stored_TMprodCost[t,dr,nd] = sum(permutedims(sum(reorder_qty,dims=1),[3,1,2]).*CpTM2); # Amount re-orderd, per supplier * the associated production cost
                # Transportation Cost TM - Supplier
                TMtransCost = reorder_delivered[:,t,:] .* CtTM';  # TM Trans Cost per time period to ship product to suppliers: Itemized out to see when costs are incurred. Rows = Suppliers Columns = Time.
                stored_TMtransCost_Itemized[:,t,dr,nd] = sum(TMtransCost,dims=2); # Tracks Transportation cost from TM to suppliers.  Itemized by supplier. Assume charge for transportation is made at the time of the order (current t), not when the items are delivered (like above variable)
                            #### CHECK IF THIS IS THE RIGHT DIMENSION
                stored_TMtransCost[t,dr,nd] = sum(sum(sum(TMtransCost,dims=2))); # Total TM trans cost for all suppliers in a single time period
                #NOTE: maybe delete on of the variables above. Or save all but last one in an auxiallary location.

                # Delivered Product 

                # Take the supplier - customer flow, reference the transportation lag for each
                # connection. Output, how much product delivered in each time period.
                # Generate a delivered matrix for each supplier that displays how much product is deliverd to each customer and in what time period.
                # determine maximum lag in the network
                temp1 = maximum(LagS);
                temp2 = maximum(LagTM);
                LagMax = maximum([temp1 temp2])
                #temp3 = maximum(LagTransferExp);
                #temp4 = maximum(LagExp);
                #LagMax = maximum([temp1 temp2 temp3 temp4]);


                ###### WORK ON
                #b = zeros(ncust,T+LagMax+LagTransferExp,G); 
                b = zeros(ncust, T+LagMax, G) # Destination matrix for delivered products. Rows = Cust Columns = Time Periods. Generate for all time periods plus largest possible lag if order was made in the last period.
                                              # Generate matrix with all customers and the total # of time periods of lag so everything combines nicely in the end

                #### NOTE: This can be done better with sums
                #Execute for Suppliers-Customers
                a =  zeros(ncust,T+LagMax,G)
                for i = 1:nsup;
                    #a = zeros(ncust,T+LagMax+LagTransferExp,G); 
                    a = zeros(ncust,T+LagMax,G)# Destination matrix same size as b for each supplier and observe when that suppliers goods are delivered to the customers.
                    for g = 1:G;
                        for j = 1:ncust;
                            a[j,t+LagS[i,j],g] = stored_transDC[i,t,g,j,dr,nd]; # 'a' is a single suppliers delivery 'windows' to all customers. 
                        end                             # Logic: Assign the flow from S1 to C1, put it in a matrix row for cust 1 and apply the time lag to when  it is delivered (column)
                    end
                    b = b .+ a; # Delivered matrix. Rows = Customers. Columns = time periods of delivery .
                end           # Combine all supplier delivery matricies into one. Final b matrix displays when the goods, that are shipped out from all suppliers in time(t), are recieved by the customers.


            """   LATERAL TRANSFER, BRING IN LATER
                  # Execute for Sup-Sup Lateral Transfers then to customers
                  # Only change is we apply a transfer lag to the transferred product
                  # only
            for i = 1:nsup;
                a = zeros(ncust,T+LagMax+LagTransferExp,G); # Destination matrix same size as b for each supplier and observe when that suppliers goods are deliverd to the customers.
                for g = 1:G;
                    for j = 1:ncust;
                        a(j,t+LagS(i,j)+LagTransferExp,g) = x.CtransferExpFLOW(i,j,g); # 'a' is a signle suppliers delivery 'windows' to all customers. 
                    end                             # Logic: Assign the flow from S1 to C1, put it in a matrix row for cust 1 and apply the time lag to when  it is delivered (column)
                end
                b = b + a ;# Delivered matrix. Rows = Customers. Columns = time periods of delivery .
            end
            """
            """ EXPEDITED SHIPPING
            # Execute for Expedited product (TM-Customer(RETAILER) direct
            for i = 1:nTM; #Should be 1x TM per product. nTM =2 would be 2x TMs per procut. KEEP AT 1
                a = zeros(ncust,T+LagMax+LagTransferExp,G); # Destination matrix same size as b for each supplier and observe when that suppliers goods are deliverd to the customers.
                for g = 1:G;
                    for j = 1:ncust;
                        a(j,t+LagExp(i,j),g) = x.CtExp(j,i,g);
                    end
                end
                b = b + a;
            end
            """
                
                stored_delivered[:,:,:,dr,nd] = stored_delivered[:,:,:,dr,nd] .+ b; # Delivered Matrix. Displays when (in terms of time) product that was shipped out in time (t) is recieved by customers. Rows = Customers Columns = Time Periods.
                            # The sum of all product in time(1).del is theamount of goods that was shipped out in time period 1.

                # Determine Backorder and Update Demand

                # Backorder
                stored_BOtotal[:,t,:,dr,nd] = dem .- stored_delivered[:,t,:,dr,nd]; # # Backorder for current time period. Rows = Customer Columns = time periods. Subtract goods delivered in current time period from demand in current time period. This is the amount of product that the customer still doesnt have, do to lag and actual shortage 
                stored_BOtrue[:,t,:,dr,nd] = value.(backorder)[:,2,:] # Rows: Customers. dimension = products. True backorder is an actual shortage of product
                stored_BOlag[:,t,:,dr,nd] =  stored_BOtotal[:,t,:,dr,nd] - stored_BOtrue[:,t,:,dr,nd] ;     
                # Want to distinguish between a time lag backorder and a shortage backorder.
                # BO is the amount of product a customer is short going from t to
                # t+1. Includes both true and time lag product
                # Time lag backorder: Inventory was on hand, it was shipped out immediatley, it just take t+x time periods for it to be delivered to customer.
                # True Shortage backorder: Not enough inventory in that time period to ship out.


                # Update demand for next period by adding the true shortage
                    if t < T # this just takes care of the edge case of the  last period. No change  to what it does.
                        dem_running[t+1,:,:] = dem_running[t+1,:,:] + stored_BOtrue[:,t,:,dr,nd];
                    else
                        #dem_running[t+1,:,:] = stored_BOtrue[:,t,:,dr,nd];
                    end
                # No need to update demand becuase backorder shortage already spoken for. It is enroute, no need to see it in the demand matrix. 
                
                
                ##### Store TC (Check on this Later)
                stored_TC[t,dr,nd] = stored_TMprodCost[t,dr,nd] + stored_TMtransCost[t,dr,nd] + stored_SupTransCost[t,dr,nd] + stored_InvCarCost[t,dr,nd] + 
                                    stored_AMprodcost[t,dr,nd] + sum(AM_force .* AM_fixed_cost)
                
            end #end of loop through t=1:60
        end #end of demand realizations
    end #end of design loop
    
    return stored_TC, stored_delivered, stored_BOtotal, stored_BOtrue, stored_BOlag, stored_AMprod, stored_transDC, stored_TMprodCost, stored_TMtransCost, stored_SupTransCost, stored_InvCarCost, stored_AMprodcost, stored_reorder_qty, stored_transDC, stored_supinv, stored_invDC, stored_reorder_delivered
end

















#################################################################################################################################################
#Only look here if you want to vary capacity by number of AM capable things. Read the next line if this applies to you.



"""   THIS ONE AND THE NEXT ONES ARE OLD. These are useful if you want to vary AM machine capacity at a single AM capable place but not over the whole system. Keeping it constant over the whole system has other benefits. This could be helpful someday. 
First MILP
This function runs the first MILP needed to solve for Big S Policy. 

Coefficient Description:
Inputs:
    T
    G
    nsup
    ncust
    CpTM
    CpAM
    CtTM
    CtS
    Ci
    CiRetailer
    CBO
    valid_TMrte
    Suprte
    CtExp
    CtransferExp
    K_stable
    LagS
    LagTM
    ndesign
    AM_fixed_Cost
    LagTransferExp
    CTransferExpFLOW
    yinit
    pAll
    forecast
    NetDesign
    AM_capacity

Used in Obj Function
    CpTM = Cost of Production Tradition Manufacturing (TM) - NumTM x T x G matrix
    CpAM = Cost of Production Additive Manufacturing (AM) - NumPossAMlocs x T x G matrix
    CtTM = Cost of Transportation from TM to Supplier (DCs) - NumTM2DCroutes x T x G matrix 
    CtS = Cost of Transportation from Supplier to SLs - NumDCs x (ncust*T) x G matrix
    Ci = Inventory Cost (at supplier) - nsup x T+1 x G matrix
    CiRetailer = Inventory Cost (at SLs) - ncsut x T+1 x G matrix
    CBO = Backorder Costs - ncust x T+1 x G matrix
Other Coefs
    Creorder - Cost of reordering (not sure if this is actually used)
    CpExp, CtExp = Expedited production + transportation costs used in Eval (TM to SL)
"""
function big_S_MILP_old(T,G,nsup,ncust,CpTM,CpAM,CtTM,CtS,Ci,CiRetailer,CBO,valid_TMrte,Suprte,CtExp,CtransferExp,K_stable,LagS,LagTM,ndesign,AM_fixed_cost,LagTransferExp,CtransferExpFLOW, yinit, pAll, forecast, NetDesign, AM_capacity)
    
    LagTM,LagS,CpTM,valid_TMrte,valid_Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable = check_dims(nsup,ncust,nTM,T,G,LagTM,LagS,CpTM,valid_TMrte,Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable)
    
    #CtExp = repeat(CtExp',1,T,G)  #Columns turn to time
    
    #CtransferExp = repeat(CtransferExp, 1, T, G) #Def needs to be modified later. Haven't tested these kinks out.
    
    ###### MIGHT BE NEEDED FOR TRANSHIPMENT BUT NEED TO WORK ON
    #valid_transferExprte = repeat(valid_transferExprte*Inf,1,T,G)
    #valid_transferExprte[isnan.(valid_transferExprte)].=0 #replace NaNs with 0
    
    #Lag Multiplier (to account for supplier lag in shipping to SLs, not a problem for DCs as current problem stands, because of one TM)
    LagSMult = ones(nsup,T,G,ncust)
    for i=1:nsup
        for j=1:ncust
            for g=1:G
                for t=1:LagS[i,j]
                    LagSMult[i,t,g,j] = 0
                end
            end
        end
    end
    K_stable = repeat(K_stable,1,T,G)

    MinPolicy = pAll #Min policy (little s)

    lbinvenDC = zeros(nsup,T+1,G)
    ubinvenDC = zeros(nsup,T+1,G) .+ (1*Inf)
    lbinvenSL = zeros(ncust,T+1,G)
    ubinvenSL = zeros(ncust, T+1,G) #.+ (Inf)
    for g=1:G
        for i=1:nsup
            lbinvenDC[i,1,g] = yinit[i,g]
            lbinvenDC[i,T+1,g] = yinit[i,g]
            ubinvenDC[i,1,g] = yinit[i,g]
            ubinvenDC[i,T+1,g] = yinit[i,g]
        end
        #### Had this loop before but it breaks it when I changed to DC perspective.
        #for i=1:ncust
            #lbinvenSL[i,1,g] = forecast[i,1,g]
            #ubinvenSL[i,1,g] = forecast[i,1,g]
        #end
    end
    
    TC_stored = zeros(ndesign)
    Big_S = zeros(nsup,G,ndesign)
    # Iterate through network designs and save the policies at the end
    for d=1:ndesign
        println("Design ", string(d))
        #println("Check Start Loop")
        AM_force = NetDesign[d,:] #Changes a network design each time (AM locations)
        K = K_stable .* AM_force # Eliminates capacity of non-AM locs. (Turn AM on/off)
        
        model1 = Model(Gurobi.Optimizer)
        #model1 = Model(GLPK.Optimizer)
        #Setting up sets for calcualtions
        NS = 1:nsup
        NC = 1:ncust
        I = 1:(nsup+ncust)
        ## Set up Decision Variables
        #println("Check DVs")
        @variable(model1, 0 <= TMprod[g=1:G,t=1:T],Int) #quantity of part g traditionally produced in period t at loc i in P
        @variable(model1, 0 <= AMprod[i=1:(nsup+ncust),t=1:T,g=1:G],Int) #quantity of part g produced with AM in period t at location i at DC or Service Loc
        @variable(model1, 0 <= transTM[ns=1:nsup,t=1:T,g=1:G] <= valid_TMrte[ns,t,g],Int) #flow of part g transported from location TM to DC in time t
        @variable(model1, 0 <= transDC[ns=1:nsup,t=1:T,g=1:G,nc=1:ncust] <= valid_Suprte[ns,t,g,nc],Int) #flow of part g from m DC loc to i SL loc at time t.
        @variable(model1, lbinvenDC[ns,t,g] <= invenDC[ns=1:nsup,t=1:T+1,g=1:G] <= ubinvenDC[ns,t,g],Int) #quantity of part g carried as inventory at location i in DC from period t-1 to t
        @variable(model1, lbinvenSL[nc,t,g] <= invenSL[nc=1:ncust,t=1:T+1,g=1:G] <= ubinvenSL[nc,t,g],Int) #quantity of part g carried as inventory at location i at SLs from period t
        @variable(model1, 0 <= backorder[nc=1:ncust,t=1:T+1,g=1:G],Int) #quantity of unsatisfied demand of part g at location i in SL from period t
        @variable(model1, 0 <= reorderQTY[ns=1:nsup,t=1:T,g=1:G],Int) #quantity to reorder for supplier i at time t for product k
        @variable(model1, 0 <= delivered[nc=1:ncust,t=1:T,g=1:G],Int) #quantity of item g "delivered/used" at SL i at time t
        @variable(model1, p[ns=1:nsup,t=1:T,g=1:G], Bin) # reordrer indicator for part g at location i (DC) at time period t
        @variable(model1, e[i=1:(ncust+nsup),t=1:T,g=1:G], Bin) #AM production indicator for part g at location i (DCs or SLs) in time period t
        @variable(model1, MinPolicy[ns,g] <= MaxPolicy[ns=1:nsup, g=1:G] <= Inf, Int) #Solving for this, Big S
        #Set up obejctive function
        @objective(model1, Min, 
            sum(CpTM[g,t]*TMprod[g,t] for g=1:G,t=1:T) + #TM production cost
            sum(CtTM[ns]*transTM[ns,t,g] for ns=1:nsup, t=1:T, g=1:G) +  #transportation cost TM to supplier
            sum(CtS[ns,nc]*transDC[ns,t,g,nc] for ns=1:nsup, t=1:T, g=1:G, nc=1:ncust) + #transportation supplier to SLs
            sum(Ci[ns,t,g]*invenDC[ns,t,g] for ns=1:nsup, t=2:T+1, g=1:G) + #DCs inventory cost
            sum(CiRetailer[nc,t,g]*invenSL[nc,t,g] for nc=1:ncust, t=2:T+1, g=1:G) + #Customer inventory cost
            sum(CpAM[i,t,g]*AMprod[i,t,g] for i=1:(nsup+ncust), t=1:T, g=1:G) + #AM production cost
            sum(CBO[nc,t,g]*backorder[nc,t,g] for nc=1:ncust, t=1:T+1, g=1:G)); #backorder cost
        
        #println("Start Constraints")
        ###Constraints
        #Flow balance produced at time T = ∑ of all things transported at time T
        #i.e. No inventory stored at TM
        for t = 1:T
            for g = 1:G
                @constraint(model1, TMprod[g,t] == sum(transTM[ns,t,g] for ns ∈ NS)) 
            end
        end

        #Inventory Balance at DCs
        for i=1:nsup
            for g=1:G
                for t=1:LagTM[i]
                    @constraint(model1, invenDC[i,t+1,g] == invenDC[i,t,g] + AMprod[i,t,g] - sum(transDC[i,t,g,nc] for nc ∈ NC))
                end
                for t=(LagTM[i]+1):T
                    @constraint(model1, invenDC[i,t+1,g] == invenDC[i,t,g] + AMprod[i,t,g] + sum(transTM[i,t-LagTM[i],g]) - sum(transDC[i,t,g,nc] for nc ∈ NC))
                end
            end
        end
        
        
        #Inventory Balance at SLs
        for i=1:ncust
            for g=1:G
                for t=1:T
                    @constraint(model1,invenSL[i,t+1,g] == invenSL[i,t,g] + AMprod[i+nsup,t,g] + 
                        sum(LagSMult[ns,t,g,i] * transDC[ns,maximum([t-LagS[ns,i],1]),g,i] for ns ∈ NS) - delivered[i,t,g])
                end
            end
        end

        #Backorder
        for i=1:ncust
            for t=1:T
                for g=1:G
                    @constraint(model1,backorder[i,t+1,g] == backorder[i,t,g] + forecast[i,t,g] - delivered[i,t,g])
                end
            end
        end

        #Min and Max Policy at DCs
        for i=1:nsup
            for g=1:G
                @constraint(model1, MinPolicy[i,g] <= invenDC[i,2,g] - (100000 * (1-p[i,1,g])))
                @constraint(model1, MaxPolicy[i,g] >= invenDC[i,2,g] - (100000 * (1-p[i,1,g])))
                for t=2:T
                    @constraint(model1, MinPolicy[i,g] <= invenDC[i,t+1,g] - (100000 * (1-p[i,t,g])) + sum(reorderQTY[i,m,g] for m=(t-1):maximum([t-LagTM[i],1])))
                    @constraint(model1, MaxPolicy[i,g] >= invenDC[i,t+1,g] - (100000 * (1-p[i,t,g])) + sum(reorderQTY[i,m,g] for m=(t-1):maximum([t-LagTM[i],1])))
                end
            end
        end

        ##### Reorder Constraints
        for i=1:nsup
            for t=1:T
                for g=1:G
                    #Reorder quantity
                    @constraint(model1, reorderQTY[i,t,g] == (MaxPolicy[i,g] - invenDC[i,t+1,g]) * p[i,t,g])

                    #Replenishment constraints
                    @constraint(model1, reorderQTY[i,t,g] <= 100000 * p[i,t,g])
                    @constraint(model1, reorderQTY[i,t,g] <= MaxPolicy[i,g] - invenDC[i,t+1,g])
                    @constraint(model1, reorderQTY[i,t,g] >= MaxPolicy[i,g] - invenDC[i,t+1,g] - 10000*(1-p[i,t,g]))
                end
            end
        end

        #Linking Reorder to trans from TM to DCs
        for i=1:nsup
            for t=1:T
                for g=1:G
                    @constraint(model1, transTM[i,t,g] == reorderQTY[i,t,g])
                end
            end
        end

        #AM capacity and production indicator
        for i=1:(nsup+ncust)
            for t=1:T
                for g=1:G
                    @constraint(model1, AMprod[i,t,g] <= K[i,t,g] * AM_capacity * e[i,t,g])
                end
                #link e with as a production indicator (wheather or not AM is currently being used for a product)
                @constraint(model1, sum(e[i,t,g] for g=1:G) <= 1)
            end
        end

        set_optimizer_attribute(model1, "OutputFlag", 1)
        set_optimizer_attribute(model1, "MIPGap", .05) #5% optimality gap
        #set_optimizer_attribute(model1, "TimeLimit", )
        set_silent(model1) #can take this out, but runs faster when you don't print the messages/results
        optimize!(model1);
        #println("Check Done")
        TC_stored[d] = objective_value(model1)

        Big_S[:,:,d] = value.(MaxPolicy)
    end
    
    return Big_S, TC_stored
end
    













""" Second MILP
    This function runs the second MILP that looks at day by day functions and makes a decision everyday then eventually calculates many Costs and returns them. 
    This is old because AM capacity here is the max capacity at a single location of AM capacity. We want to keep this constant for research purposes. 
"""
function MILP2_old(T,G,nsup,ncust,CpTM,CpAM,CtTM,CtS,Ci,CiRetailer,CBO,valid_TMrte,Suprte,CtExp,CtransferExp,K_stable,LagS,LagTM,ndesign,AM_fixed_cost,LagTransferExp,CtransferExpFLOW, yinit, pAll, dem_real, NetDesign, AM_capacity, Big_S) 
    #### Eval Code
    #nd = 16  #Use this if testing single runs to check for speed 
    #initialize all stored variables for later
    maxT = maximum([maximum(LagS) maximum(LagTM)]) #some arrays need to be past time horizon because of shipping lag, this is the maxT to save space later
    stored_transDC = zeros(nsup,T,G,ncust,nreal,ndesign);
    stored_AMprod = zeros(nsup+ncust,T,G,nreal,ndesign);
    stored_SupTransCost = zeros(T,nreal,ndesign);
    stored_AMprodcost = zeros(T,nreal,ndesign);
    stored_TMprodCost = zeros(T,nreal,ndesign);
    stored_TMprod = zeros(T,nreal,ndesign);
    stored_InvCarCost = zeros(T,nreal,ndesign);
    stored_TC = zeros(T,nreal,ndesign);
    stored_reorder_qty = zeros(nsup,T,G,nreal,ndesign);
    stored_TMtransCost_Itemized = zeros(nsup,T,nreal,ndesign);
    stored_delivered = zeros(ncust,T+maxT,G,nreal,ndesign);  #####may need to change the max later
    #stored_reorder = zeros();
    stored_reorder_delivered = zeros(nsup,T+maxT,G,nreal,ndesign); 
    stored_TMtransCost = zeros(T,nreal,ndesign);
    stored_BOtotal = zeros(ncust,T+1,G,nreal,ndesign);
    stored_BOtrue = zeros(ncust, T+1, G, nreal,ndesign);
    stored_BOlag = zeros(ncust,T+1,G,nreal,ndesign);
    stored_supinv = zeros(nsup,T+7,G,nreal,ndesign);
    stored_invDC = zeros(nsup,T+7,G,nreal,ndesign);


    LagTM,LagS,CpTM,valid_TMrte,valid_Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable = check_dims(nsup,ncust,nTM,T,G,LagTM,LagS,CpTM,valid_TMrte,Suprte,Ci,CiRetailer,CBO,CpAM,AM_fixed_cost,K_stable)
    
    #### if not solving take "@distributed" out and error should pop up. 
    for nd = 1:ndesign
        #update AM force
        println("Design ", string(nd))
        AM_force = NetDesign[nd,:] #changes the network design of AM lcoations

        #Setting Min/Max Policy found in MILP
        MaxPolicy = Big_S[:,:,nd]
        if nsup == 1
            MaxPolicy_3d = reshape(MaxPolicy,(nsup,1,G)) #changes to be rows=nsup, cols=t, 3d=product
            MinPolicy = pAll
            MinPolicy_3d = reshape(MinPolicy,(nsup,1,G))
        else
            MaxPolicy_3d = reshape(MaxPolicy,(nsup,1,G)) #changes to be rows=nsup, cols=t, 3d=product
            MinPolicy = pAll
            MinPolicy_3d = reshape(MinPolicy,(nsup,1,G))
        end
            

        #Loop for demand realizations
        for dr = 1:nreal
            dem_running = dem_real[dr,:,:,:]
            #needs to be reset every run 
            #starting with 5 more products than little s
            yinit = MinPolicy .+ 5
            
            #setting up inventories
            sup = zeros(nsup,T+maxT,G)
            #setup so sup is nsup, time, by nprods matrix,
            for i = 1:nsup
                for g = 1:G
                    sup[i,1,g] = yinit[i,g]
                    stored_supinv[i,1,g,dr,nd] = yinit[i,g]
                end
            end


            previous_orders = zeros(nsup,T+maxT,G) # reset for each realization
            reorder_delivered = zeros(nsup,T+maxT,G) # Destination  matrix for when the re-ordered product will be delivered. Rows = Suppliers. Columns = Time

            #iterate through time
            for t=1:T
                #println("This is time period ", string(t))
                ## Reduce dimensions of select variables for evaluation in MILP
                # Assumption CtS, CBO, CpAM, K do not change over time. If they do then
                # just need to extract the appropriate columns on the appropriate loop
                # through time.
                dem = dem_running[t,:,:];# extract the current time period for use in MILP
                Ci2 = Ci[:,t:t+1,:]; # reduce to 2
                Ci2[:,1,:] .= 0 # First time period (yinit = 0 cost)
                CBO2 = CBO[:,t:t+1,:]; # reduce to 2
                CBO2[:,1,:] .= 0;

                CpAM2 = CpAM[:,t,:] # reduce to 1

                # Reduce to 1
                # Eliminate Capacity of Non-AM 'forced' Sites
                K = K_stable.*AM_force;


                CiRetailer2 = CiRetailer[:,t:t+1,:]; # reduce to 2
                CiRetailer2[:,1,:] .= 0;



                # Not used in MILP but used in future calcs
                CpTM2 = CpTM[:,t,:];

                # Reduce Expedited TM-Cust
                #CpExp2 = CpExp[:,t,:]; # reduce to 1
                #CtExp2 = CtExp[:,t,:]; # reduce to 1

                #CtExpCapacity2 = CtExpCapacity[:,t,:];


                # Reduce the valid rte UBs for use in MILP
                valid_Suprte2 = valid_Suprte[:,t,:,:]; # Used in MILP
                valid_TMrte2 = valid_TMrte[:,t,:]; # Not used in this code..


                # Reduce Sup-Sup Exp transfer variables and UBs
                #CtransferExp2 = CtransferExp[:,((t-1)*3+1):(t-1)*3+3,:];
                #valid_transferExprte2 = valid_transferExprte[:,((t-1)*3+1):(t-1)*3+3,:];

                #CtransferExpFLOW2 = CtransferExpFLOW[:,t:t+1,:];

                yinit = sup[:,t,:] # update inventory for current t
                lbSupply = zeros(nsup,2,G)
                ubSupply = ones(nsup,2,G) .* Inf
                for i = 1:nsup
                    for g=1:G
                        lbSupply[i,1,g] = yinit[i,g]
                        ubSupply[i,1,g] = yinit[i,g]
                    end
                end
                #yinit = permute(yinit,[2 1 3]); 
                y = zeros(nsup,T,G)

                #Set up sets (not necessarily needed but is used)
                NC = 1:ncust
                NS = 1:nsup
                
                #println("CHECK FOR YINIT AND LBSUP AND UBSUP")
                #println(yinit)
                #println(lbSupply)
                #println(ubSupply)
                
                #model = Model(Gurobi.Optimizer);
                model = Model(GLPK.Optimizer)
                ## Set up Decision Variables
                 #mp = buildMILP(Ci,CpAM,CtS,cdel,CiRetailer,CBO,e,CpExp,CtExp,CtransferExp,CtransferExpFLOW)
                @variable(model, lbSupply[ns,tt,g] <= invenDC[ns=1:nsup,tt=1:2,g=1:G] <= ubSupply[ns,tt,g]) #inventory quantity of prod g at DC nc at time t
                @variable(model, 0 <= AMprod[i=1:(nsup+ncust),g=1:G] <= K[i]*AM_capacity) #quantity of product g produced at location i
                #@variable(model, 0 <= transTM[i=1:nsup,g=1:G]) 
                @variable(model, 0 <= transDC[ns=1:nsup,g=1:G,nc=1:ncust] <= valid_Suprte2[ns,g,nc]) #quantity of product g transported from DC ns to SL nc
                @variable(model, 0 <= delivered[nc=1:ncust, g=1:G]) #Quantity of product g delivered/utilized at SL nc
                @variable(model, 0 <= invenSL[nc=1:ncust,tt=1:2,g=1:G] <= 0) #Inventory at SL
                @variable(model, 0 <= backorder[nc=1:ncust,tt=1:2,g=1:G]) #Backorder DV 
                @variable(model, e[i=1:(ncust+nsup),g=1:G], Bin) #AM production indicator
                #@variable(model, 0 <= ExpProd[i=1:ncust,g=1:G]) #Expedited Production
                #@variable(model, 0 <= transExp[i=1:ncust,g=1:G]) #Expedited transportation
                #variable(model, transferExp[]) #Lateral Transfer between DCs

                #Set up obejctive function
                @objective(model, Min, 
                    sum(Ci2[ns,tt,g] * invenDC[ns,tt,g] for ns=1:nsup, tt=1:2, g=1:G) + #DC inventory cost
                    sum(CpAM2[i,g] * AMprod[i,g] for i=1:(nsup+ncust), g=1:G) + #AM production cost
                    sum(CtS[ns,nc] * transDC[ns,g,nc] for ns=1:nsup, g=1:G, nc=1:ncust) + #DC to SL transportation cost
                    sum(CiRetailer2[nc,tt,g] * invenSL[nc,tt,g] for nc=1:ncust, tt=1:2, g=1:G) + #SL inventory cost
                    sum(CBO2[nc,tt,g] * backorder[nc,tt,g] for nc=1:ncust, tt=1:2, g=1:G)) #+ #Backorder Cost
                    #sum(CpExp2[nc,g] * ExpProd[nc,g] for nc=1:ncust, g=1:G) + #Expedited production cost
                    #sum(CtransferExp2[nc,g] * transExp[nc,g] for nc=1:ncust, g=1:G) + #Expedtied transportation cost


                for g = 1:G
                    for tt=1:1
                        # Expedited TM Flow Balance
                        # Amount Exp Produced - sum of Amount EXPEDITED out to all RETAILERS.
                        #@constraint(model, ExpProd[g] == sum(transEXP[i,g] for i=1:nsup))

                        for ns = 1:nsup

                            # Distribution Center (Supplier) Flow Balance
                            # In from TM, + Current Inv + AM Produced - Leftover Inv - Sum of everything shipped out. 
                            @constraint(model, invenDC[ns,tt,g] + AMprod[ns,g] - invenDC[ns,tt+1,g] - sum(transDC[ns,g,nc] for nc ∈ NC) == 0 )
                            # *** want to add: in all in 1st column, all rows and all out: 1st row all columns


                            # Build the transfer connections for 'In and out'. It is OK if the diagnoal is 'wrong' bc the UB on it is zero so nothing will ever flow over it.       
                            #zz = zeros(size(CtransferExp));
                            #zz[:,ns,g] = 1;
                            #zz[ns,:,g] = -1;
                            #mp.addcstr({[1 -1],{ns,[tt tt+1],g}},{ns,tt,g},{-1,{ns,':',g}},0,0,0,0,0,0,{zz,{':',':',g}},{-1,{ns,':',g}},'=',0);

                            #mp.addcstr(0,0,0,0,0,0,0,0,0,{':',ns,g},'=',{ns,':',g}); # Everything in from a lateral transfer must be transported across the special transfer arcs.


                            end #end ns

                        for nc = 1:ncust;
                            # Retailer Flow Balance
                            # Sum of everything in from all suppliers - delivered +Inv In - Inv out
                                # This places a '1' in the appropraite time period that product left the supplier based on the lag.
                            @constraint(model, AMprod[nc+nsup,g] + sum(transDC[ns,g,nc] for ns ∈ NS) - delivered[nc,g] + invenSL[nc,tt,g] - invenSL[nc,tt+1,g] == 0)

                            #BO
                            # in from delivered, -BO in + BO out, = Demand
                            @constraint(model, delivered[nc,g] - backorder[nc,tt,g] + backorder[nc,tt+1,g]== dem[nc, g])
                            end #end nc

                        # Production indicator constraint and AM capacity enforecement
                        #  Fixed Cost no longer active in MILP.
                        # ****************indedxing is different becuase it is for both
                        # suppliers and customers.**********
                        for i = 1:nsup+ncust;
                            @constraint(model, AMprod[i,g] <= K[i] * e[i,g] * AM_capacity)
                        end
                    end #end of tt loop
                end #end of G loop
                #link e with as a production indicator (wheather or not AM is currently being used for a product)
                for i = 1:(nsup+ncust)
                    @constraint(model, sum(e[i,g] for g=1:G) <= 1)
                end
                set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_ALL)
                #set_optimizer_attribute(model, "OutputFlag", 1)
                #set_optimizer_attribute(model, "MIPGap", .001)
                #set_optimizer_attribute(model, "TimeLimit", )
                set_silent(model) #can take this out, but runs faster when you don't print the messages/results
                optimize!(model);

                #store the data 
                stored_transDC[:,t,:,:,dr,nd] .= value.(transDC)
                stored_AMprod[:,t,:,dr,nd] .= value.(AMprod)
                stored_SupTransCost[t,dr,nd] = sum(CtS[ns,nc] * value.(transDC)[ns,g,nc] for ns=1:nsup, g=1:G, nc=1:ncust)
                stored_AMprodcost[t,dr,nd] = sum(CpAM2[i,g] * value.(AMprod)[i,g] for i=1:(nsup+ncust), g=1:G)
                stored_InvCarCost[t,dr,nd] = sum(Ci2[ns,tt,g] * value.(invenDC)[ns,tt,g] for ns=1:nsup, tt=1:2, g=1:G)
                stored_invDC[:,t+1,:,dr,nd] = value.(invenDC)[:,2,:]
                
                #println("Value of invenDC")
                #println(value.(invenDC))
                #println(size(value.(invenDC)))
                

                
                #Update inventory for next time period
                sup[:,t+1,:] = value.(invenDC)[:,2,:] + sup[:,t+1,:] # Represents the amount of inventory on hand to START that period.  Row = Suppliers. Columns = Time. 
                stored_supinv[:,t+1,:,dr,nd] = sup[:,t+1,:]
                ######MIGHT BE OF INTEREST LATER NEED TO BUILD STRUCTURE FOR THIS STORAGE
                ## AM locations used
                #AM_loc_used(:,t,:) = x.CpAM > 0;
                ## AM locations used
                #AM_qty_used(:,t,:) = x.CpAM;

                #Inventory Carried
                
                
                y[:,t,:] = value.(invenDC)[:,2,:] # Inventory carried between periods. (t1 is the inventory that was carried from t1 to t2)
                                # This is different than 'sup' at the end of the loop becuase y extracts your inventory before re-orders arrive.
                                # Use y for inventory carried between periods.
                                # Use sup to see what the actual inventory is at the start of a period (it adds in re-orders)

                # Determine the re-order quantity and place the re-ordered quantiy
                # in sup matrix based on the transportation lag from TM to sup.

                # This should prevent re-order every period, by looking back over
                # the lag time. If your inventory + whatever you have enroute is <
                # Minpolicy then order again
                
                log = zeros(nsup,1,G)
                if t==1
                    log = reshape((y[:,t,:] .+ previous_orders[:,t,:] .< MinPolicy),(nsup,1,G))
                    #println("Inv and ordered")
                    #println((y[:,t,:] .+ previous_orders[:,t,:]))
                else 
                    log = reshape((y[:,t,:] .+ permutedims(sum(previous_orders[:,maximum([1, t-6]):t-1,:],dims=2),[1,3,2])) .< MinPolicy,(nsup,1,G))
                    #println("Inv and order")
                    #println(y[:,t,:] .+ permutedims(sum(previous_orders[:,maximum([1, t-6]):t-1,:],dims=2),[1,3,2]))
                end
                ################### if you change LagTM to non-uniform look at if statement above and edit accordingly (change the 6 somehow), most likely need a loop for each product (if change lag)
 
                
                b = log .* MaxPolicy_3d; # Zeros out corresponding policies that dont need to order
                bb = log .* reshape(y[:,t,:],(nsup,1,G)); # Zeros out corresponding suppliers that dont need to order
                reorder_qty  = b.-bb; # subtract current inventory from policies. This is the re-order qty per supplier.
                previous_orders[:,t,:] = reorder_qty
                
                
                # Now place in the appropriate place (using nested loop to get it right) ##COMMENT: g loop might not be needed because same lag time
                for ns=1:nsup
                    for g=1:G
                        sup[ns,t+LagTM[ns],g] = reorder_qty[ns,1,g] + sup[ns,t+LagTM[ns],g] # Add the re-ordered qty to the location where it will be delivered.
                        stored_supinv[ns,t+LagTM[ns],g,dr,nd] = sup[ns,t+LagTM[ns],g]
                        reorder_delivered[ns,t+LagTM[ns],:] = reorder_qty[ns,1,:] # Just a data tracking matrix. Displays when re-ordered goods, for  each supplier  are  delivered.
                        stored_reorder_delivered[ns,t+LagTM[ns],:,dr,nd] = reorder_qty[ns,1,:]
                    end
                end

                # Save Data  
                stored_reorder_qty[:,t,:,dr,nd] = reorder_qty ; # Reordered quantity per customer per product
                

                """ NOT SURE IF NEEDED

                # NOTE: Sup matrix recieved the re-ordered goods.  At this point,     ################ WORK ON
                # it represents the full amount of inventory available at
                # suppliers for t+1 and on.
                time(t).sup = sup; # This is the inventory available at the beginiing of the next time period time period.
                              # ex. time(1).sup represents the supply conditions to start t2.
                              # NOTE: AM Capacity is added in at the beginning
                              # of next period.
                """

                #Costs associated with re-orders
                
                # Production Cost TM
                stored_TMprodCost[t,dr,nd] = sum(permutedims(sum(reorder_qty,dims=1),[3,1,2]).*CpTM2); # Amount re-orderd, per supplier * the associated production cost
                # Transportation Cost TM - Supplier
                TMtransCost = reorder_delivered[:,t,:] .* CtTM';  # TM Trans Cost per time period to ship product to suppliers: Itemized out to see when costs are incurred. Rows = Suppliers Columns = Time.
                stored_TMtransCost_Itemized[:,t,dr,nd] = sum(TMtransCost,dims=2); # Tracks Transportation cost from TM to suppliers.  Itemized by supplier. Assume charge for transportation is made at the time of the order (current t), not when the items are delivered (like above variable)
                            #### CHECK IF THIS IS THE RIGHT DIMENSION
                stored_TMtransCost[t,dr,nd] = sum(sum(sum(TMtransCost,dims=2))); # Total TM trans cost for all suppliers in a single time period
                #NOTE: maybe delete on of the variables above. Or save all but last one in an auxiallary location.

                # Delivered Product 

                # Take the supplier - customer flow, reference the transportation lag for each
                # connection. Output, how much product delivered in each time period.
                # Generate a delivered matrix for each supplier that displays how much product is deliverd to each customer and in what time period.
                # determine maximum lag in the network
                temp1 = maximum(LagS);
                temp2 = maximum(LagTM);
                LagMax = maximum([temp1 temp2])
                #temp3 = maximum(LagTransferExp);
                #temp4 = maximum(LagExp);
                #LagMax = maximum([temp1 temp2 temp3 temp4]);


                ###### WORK ON
                #b = zeros(ncust,T+LagMax+LagTransferExp,G); 
                b = zeros(ncust, T+LagMax, G) # Destination matrix for delivered products. Rows = Cust Columns = Time Periods. Generate for all time periods plus largest possible lag if order was made in the last period.
                                              # Generate matrix with all customers and the total # of time periods of lag so everything combines nicely in the end

                #### NOTE: This can be done better with sums
                #Execute for Suppliers-Customers
                a =  zeros(ncust,T+LagMax,G)
                for i = 1:nsup;
                    #a = zeros(ncust,T+LagMax+LagTransferExp,G); 
                    a = zeros(ncust,T+LagMax,G)# Destination matrix same size as b for each supplier and observe when that suppliers goods are delivered to the customers.
                    for g = 1:G;
                        for j = 1:ncust;
                            a[j,t+LagS[i,j],g] = stored_transDC[i,t,g,j,dr,nd]; # 'a' is a single suppliers delivery 'windows' to all customers. 
                        end                             # Logic: Assign the flow from S1 to C1, put it in a matrix row for cust 1 and apply the time lag to when  it is delivered (column)
                    end
                    b = b .+ a; # Delivered matrix. Rows = Customers. Columns = time periods of delivery .
                end           # Combine all supplier delivery matricies into one. Final b matrix displays when the goods, that are shipped out from all suppliers in time(t), are recieved by the customers.


            """   LATERAL TRANSFER, BRING IN LATER
                  # Execute for Sup-Sup Lateral Transfers then to customers
                  # Only change is we apply a transfer lag to the transferred product
                  # only
            for i = 1:nsup;
                a = zeros(ncust,T+LagMax+LagTransferExp,G); # Destination matrix same size as b for each supplier and observe when that suppliers goods are deliverd to the customers.
                for g = 1:G;
                    for j = 1:ncust;
                        a(j,t+LagS(i,j)+LagTransferExp,g) = x.CtransferExpFLOW(i,j,g); # 'a' is a signle suppliers delivery 'windows' to all customers. 
                    end                             # Logic: Assign the flow from S1 to C1, put it in a matrix row for cust 1 and apply the time lag to when  it is delivered (column)
                end
                b = b + a ;# Delivered matrix. Rows = Customers. Columns = time periods of delivery .
            end
            """
            """ EXPEDITED SHIPPING
            # Execute for Expedited product (TM-Customer(RETAILER) direct
            for i = 1:nTM; #Should be 1x TM per product. nTM =2 would be 2x TMs per procut. KEEP AT 1
                a = zeros(ncust,T+LagMax+LagTransferExp,G); # Destination matrix same size as b for each supplier and observe when that suppliers goods are deliverd to the customers.
                for g = 1:G;
                    for j = 1:ncust;
                        a(j,t+LagExp(i,j),g) = x.CtExp(j,i,g);
                    end
                end
                b = b + a;
            end
            """
                
                stored_delivered[:,:,:,dr,nd] = stored_delivered[:,:,:,dr,nd] .+ b; # Delivered Matrix. Displays when (in terms of time) product that was shipped out in time (t) is recieved by customers. Rows = Customers Columns = Time Periods.
                            # The sum of all product in time(1).del is theamount of goods that was shipped out in time period 1.

                # Determine Backorder and Update Demand

                # Backorder
                stored_BOtotal[:,t,:,dr,nd] = dem .- stored_delivered[:,t,:,dr,nd]; # # Backorder for current time period. Rows = Customer Columns = time periods. Subtract goods delivered in current time period from demand in current time period. This is the amount of product that the customer still doesnt have, do to lag and actual shortage 
                stored_BOtrue[:,t,:,dr,nd] = value.(backorder)[:,2,:] # Rows: Customers. dimension = products. True backorder is an actual shortage of product
                stored_BOlag[:,t,:,dr,nd] =  stored_BOtotal[:,t,:,dr,nd] - stored_BOtrue[:,t,:,dr,nd] ;     
                # Want to distinguish between a time lag backorder and a shortage backorder.
                # BO is the amount of product a customer is short going from t to
                # t+1. Includes both true and time lag product
                # Time lag backorder: Inventory was on hand, it was shipped out immediatley, it just take t+x time periods for it to be delivered to customer.
                # True Shortage backorder: Not enough inventory in that time period to ship out.


                # Update demand for next period by adding the true shortage
                    if t < T # this just takes care of the edge case of the  last period. No change  to what it does.
                        dem_running[t+1,:,:] = dem_running[t+1,:,:] + stored_BOtrue[:,t,:,dr,nd];
                    else
                        #dem_running[t+1,:,:] = stored_BOtrue[:,t,:,dr,nd];
                    end
                # No need to update demand becuase backorder shortage already spoken for. It is enroute, no need to see it in the demand matrix. 
                
                
                ##### Store TC (Check on this Later)
                stored_TC[t,dr,nd] = stored_TMprodCost[t,dr,nd] + stored_TMtransCost[t,dr,nd] + stored_SupTransCost[t,dr,nd] + stored_InvCarCost[t,dr,nd] + 
                                    stored_AMprodcost[t,dr,nd] + sum(AM_force .* AM_fixed_cost)
                
            end #end of loop through t=1:60
        end #end of demand realizations
    end #end of design loop
    
    return stored_TC, stored_delivered, stored_BOtotal, stored_BOtrue, stored_BOlag, stored_AMprod, stored_transDC, stored_TMprodCost, stored_TMtransCost, stored_SupTransCost, stored_InvCarCost, stored_AMprodcost, stored_reorder_qty, stored_transDC, stored_supinv, stored_invDC, stored_reorder_delivered
end