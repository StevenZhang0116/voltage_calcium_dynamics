NEURON {
        SUFFIX steady_k
        USEION k WRITE ik VALENCE 1
        RANGE ik, rate
}

PARAMETER {
    rate    (mA/cm2)
}

INITIAL {
    rate = 1
}
 
ASSIGNED {
    ik     (mA/cm2)
}
 
BREAKPOINT {
    ik = rate
}

