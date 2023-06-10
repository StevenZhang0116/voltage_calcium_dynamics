function outFreq = compSig(nSoma, nDend)
    a1 = 87.01;
    a2 = 68.24;
    a3 = 71.71;
    a4 = 10.5;
    b1 = 28.5;
    b2 = 164.7;
    b3 = 64.97;
    b4 = -12.63;
    M = a1 + a2./(1+exp(-(nDend-a3)./a4));
    T = b1 + b2./(1+exp(-(nDend-b3)./b4));
    outFreq = M./(1+exp(-(nSoma-T)));
end