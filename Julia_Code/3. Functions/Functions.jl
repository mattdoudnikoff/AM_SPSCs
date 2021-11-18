"""
Function takes a transition matrix and specified nodes of the matrix and returns number of each kind of node as 
well as path valid paths from TM to DC and DC to SL.

Inputs: 
    trans_matrix: total node X total node matrix with valid paths. (Contains a 1 or nonzero number for valid path. 0 o/w)
    TM: vector of which nodes are TMs
    DC: vector of which nodes are DCs
    SL: vector of which nodes are SLs
Outputs:
    nsup: number of suppliers
    ncust: number of customers
    nTM: number of TMs
    valid_TMrte: nsup X nTM X G matrix that is 1 for valid route or 0 otherwise for each product from TM to DCs #####NOTE THIS PROBABLY WILL CHANGE LATER
    Suprte: nsup X ncust marix that is a 1 for valid route and 0 otherwise.
"""
function trans_matrix_math(trans_matrix, TM, DC, SL)
    #Get the number of each kind of node
    nsup = length(DC)
    ncust = length(SL)
    nTM = length(TM)

    #Build a matrix that contains that routes and costs between TM to DC (is a nsup X nTM X G matrix)
    valid_TMrte = zeros(nsup,nTM,G)
    CtTM = zeros(nTM, nsup)
    for i=1:nTM
        valid_TMrte[:,i,:] .= (trans_matrix[TM[i],DC] .> 0)'
        CtTM[i,:] = trans_matrix[TM[i],DC]
    end
    
    #Biuld a matrix that has routes between DC to SL (nsup x ncust matrix)
    Suprte = zeros(nsup,ncust)
    CtS = zeros(nsup,ncust)
    for i=1:nsup
        Suprte[i,:] = (trans_matrix[DC[i],SL] .> 0)
        CtS[i,:] = trans_matrix[DC[i],SL]
    end
    
    
    
    return nsup, ncust, nTM, valid_TMrte, Suprte, CtTM, CtS
end

"""
Function to calcualte 95% Confidence Interval based on list 1d array of data.
"""
function CI95(data)
    z = 1.96
    lower = mean(data) - (z*std(data)/sqrt(length(data)))
    upper = mean(data) + (z*std(data)/sqrt(length(data)))
    return lower, upper
end

"""
Function that takes the mess of many dimensional outputs from MILP2/Eval phase and gives them back in a prettier/easier way to see it (needs work but decent for now)
"""
function analyze_outputs(TC, delivered, BOtot, AMprod, transDC, TMprodCost, TMtransCost, SupTransCost, InvCarCost, AMprodCost, demand)
    avg_TC_per_Des = permutedims(mean(sum(TC,dims=1),dims=2),[3,1,2])
    TCstdev = permutedims(std(sum(TC,dims=1),dims=2),[3,1,2])
    
    avg_TMprodCost_per_Des = permutedims(mean(sum(TMprodCost,dims=1),dims=2),[3,1,2])
    TMcoststdev = permutedims(std(sum(TMprodCost,dims=1),dims=2),[3,1,2])
    
    avg_TMtransCost_per_Des = permutedims(mean(sum(TMtransCost,dims=1),dims=2),[3,1,2])
    TMtransstdev = permutedims(std(sum(TMtransCost,dims=1),dims=2),[3,1,2])
    
    avg_SupTransCost_per_Des = permutedims(mean(sum(SupTransCost,dims=1),dims=2),[3,1,2])
    Suptransstdev = permutedims(std(sum(SupTransCost,dims=1),dims=2),[3,1,2])
    
    avg_InvCarCost_per_Des = permutedims(mean(sum(InvCarCost,dims=1),dims=2),[3,1,2])
    InvCoststdev = permutedims(std(sum(InvCarCost,dims=1),dims=2),[3,1,2])
    
    avg_AMprodCost_per_Des = permutedims(mean(sum(AMprodCost,dims=1),dims=2),[3,1,2])
    AMprodcoststdev = permutedims(std(sum(AMprodCost,dims=1),dims=2),[3,1,2])
    
    avg_BO_by_SL = permutedims(mean(sum(sum(BOtot,dims=2),dims=3),dims=4),[5,1,2,3,4])
    avg_BO_per_Des = permutedims(mean(sum(sum(sum(BOtot,dims=1),dims=2),dims=3),dims=4),[5 1 2 3 4])
    BOstdev = permutedims(std(sum(sum(sum(BOtot,dims=1),dims=2),dims=3),dims=4),[5 1 2 3 4])
    avg_AMprod_per_Des = mean(sum(sum(sum(AMprod,dims=1),dims=2),dims=3),dims=4)
    percentByAM = permutedims(mean(permutedims(sum(sum(sum(AMprod,dims=1),dims=2),dims=3),[4,5,1,2,3]) ./ sum(sum(sum(demand,dims=2),dims=3),dims=4),dims=1),[2,1,3,4,5])
    return avg_TC_per_Des, avg_BO_per_Des, avg_AMprod_per_Des, avg_TMprodCost_per_Des, 
        avg_TMtransCost_per_Des, avg_SupTransCost_per_Des,avg_InvCarCost_per_Des, avg_AMprodCost_per_Des, percentByAM, TCstdev, BOstdev, avg_BO_by_SL, TMcoststdev, TMtransstdev, Suptransstdev, InvCoststdev, AMprodcoststdev
end

"""
Function to save data appropriately in the right spot with a description
    NOTE: THIS IS ALWAYS A WORK IN PROGRESS AND CHANGE AS YOU NEED TO SAVE OTHER VARIABLES AND COSTS AND THINGS OF INTEREST
"""
function save_data(description,adi,cv2,avgdem,TC, BOtot, AMproduced, TMprodCost, TMtransCost, SupTransCost, InvCarCost, AMprodCost, Big_S, little_s, percentAM, TCstdev, BOstdev, BO_by_SL, TMcoststdev, TMtransstdev, Suptransstdev, InvCoststdev, AMprodcoststdev)
    thispath = pwd()
    if thispath[1] == '/'
        added = "9. Raw Output Data/"
    else
        added = "9. Raw Output Data\\"
    end
    date = string(today())
    adi_str = replace(string(adi),"." => "_")
    cv2_str = replace(string(cv2),"." => "_")
    avgdem_str = string(avgdem)
    path = thispath[1:findlast("Julia_Code",thispath)[end]+1] * added * date * "_" * description * "_" * adi_str * "ADI_" * cv2_str * "CV2" * "_" * avgdem_str * "AvgDem" * ".jld2"
    println(path)
    @save path TC BOtot AMproduced TMprodCost TMtransCost SupTransCost InvCarCost AMprodCost Big_S little_s percentAM TCstdev BOstdev BO_by_SL TMcoststdev TMtransstdev Suptransstdev InvCoststdev AMprodcoststdev
    return "Data saved here: " * path
end


