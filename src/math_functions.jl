function welford_online(M, S, n, x)
    oldM = M
    M = M + (x - M) / n
    S = S + (x - M) * (x - oldM)
    return M, S 
end