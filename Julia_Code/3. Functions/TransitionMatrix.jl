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