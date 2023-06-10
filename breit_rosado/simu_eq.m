function dydt = simu_eq(y,impulse)
    % y = [c1,o2,c2,o1,c_c,b]

    % define constants
    ka_minus = 28.8; % s-1
    ka_plus = 1500; % \mu M-4 s-1
    kb_minus = 385.9; % s-1
    kb_plus = 1500; % \mu M-3 s-1
    kc_minus = 0.1; % s-1
    kc_plus = 1.75; % s-1
    
    btot = 40; % \mu M
    D_c = 220; % \mu m2 s-1
    D_b = 20; % \mu m2 s-1
    kapb_plus = 27; % \mu M-1 s-1
    kapb_minus = 19; % s-1

    % given impulse
    y(5) = y(5)+impulse; 
    
    % solve eq
    dydt = zeros(4,1);
    % c_1 Eq 15
    dydt(1) = ka_minus*y(4)-ka_plus*y(5)^4*y(1);
    % o_2 Eq 16
    dydt(2) = kb_plus*y(5)^3*y(4)-kb_minus*y(2);
    % c_2 Eq 17
    dydt(3) = kc_plus*y(4)-kc_minus*y(3);
    % o_1 Eq 14
    dydt(4) = -(dydt(1)+dydt(2)+dydt(3));
    % c_c Eq 3, Laplacian=0
    dydt(5) = (kapb_minus*(btot-y(6))-kapb_plus*y(6)*y(5));
    % b Eq 4, Laplacian=0
    dydt(6) = (kapb_minus*(btot-y(6))-kapb_plus*y(6)*y(5));
   
end