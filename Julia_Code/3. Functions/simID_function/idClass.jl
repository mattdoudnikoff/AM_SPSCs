using Statistics
function idClass(data)
    #data: time series demand data. Rows: time series. Columns: time periods
    # p: average inter-demand interval (per time series)
    # cv2: squared coeficient of variance (per time series)
    # lev: mean level of demand (per time series)
    # ***_avg: average value of all time series combined
    N = size(data,1);
    
    #initiate variables
    p = zeros(N)
    cv2 = zeros(N)
    lev = zeros(N)
    for n in 1:N
        ts = data[n,:]; # extracts a single time series
        nzd = findall(ts .!= 0); #idx of non-zero demand
        k = length(nzd); # number of non-zero demand periods
        z = ts[nzd]; # non-zero demand values
        x = hcat(nzd[1], (nzd[2:k] - nzd[1:(k-1)])'); # inter-demand interval
        p[n] = mean(x); # adi
        cv2[n] = (std(z)/mean(z))^2; # cv2
        lev[n] = mean(z); # lev
    end
    # Consolidate all time series
    p_avg = mean(p);
    cv2_avg = mean(cv2);
    lev_avg = mean(lev);

    # Classify forecasting method
    use_SBA = ((p .> 1.32) .| (cv2 .> 0.49)); # use SBA for these scenarios
    use_Croston = (use_SBA .== 0); # use croston for everything else (bottom left corner)

    # Find idx of demand streams that go with each forecast type
    idx_SBA = findall(use_SBA);
    idx_Croston = findall(use_Croston);
    
    return [p, cv2, lev, p_avg, cv2_avg, lev_avg, idx_SBA, idx_Croston]
end