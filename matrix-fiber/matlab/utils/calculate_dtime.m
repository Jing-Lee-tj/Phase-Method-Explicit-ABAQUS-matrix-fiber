function dtime = calculate_dtime(E1,v12,rho,L_min,c0,eta_L,l0_L,eta_m,l0_m)
    s = 0.6; %security factor
    lamda_ = E1 * v12 / ( (1+v12) * (1-2*v12));
    miu_ = E1 / ( 2* (1+v12));
    c_d = sqrt( (lamda_ + 2*miu_) / rho );
    dtime_u = calculate_dtime_u(L_min,c_d);
    dtime_L = calculate_dtime_L(L_min,c0,eta_L,l0_L);
    dtime_m = calculate_dtime_m(L_min,c0,eta_m,l0_m);
    dtime = s * min([dtime_u,dtime_L,dtime_m]);
end

function dtime = calculate_dtime_u(L_min,c_d)
    dtime = L_min/c_d;
end
function dtime = calculate_dtime_L(L_min,c0,eta_L,l0_L)
    dtime = c0 * eta_L * L_min * L_min /4 / l0_L;
end
function dtime = calculate_dtime_m(L_min,c0,eta_m,l0_m)
    beita = 35; %penalty parameter
    dtime = c0 * eta_m * L_min * L_min /4 / l0_m / (1+beita);
end