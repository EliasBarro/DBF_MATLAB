% Michael Qiu - 2022

function V_max = max_cruiseV(T_max,Cd_min,AR,W,S,ALT)

    k = liftind_CD(AR);
    rho = dens(ALT);
    
    V_max = sqrt((T_max+sqrt(T_max.^2-4*Cd_min*k.*W.^2))/(rho.*S.*Cd_min)); %% See Snorri pg. 877, Eqn. 19-25

end