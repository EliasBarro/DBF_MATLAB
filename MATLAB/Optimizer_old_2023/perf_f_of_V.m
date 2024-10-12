% Michael Qiu - 2022

function f = perf_f_of_V(ALT,S,Cd_min,AR,W,P,eff,V)

rho = dens(ALT);

k = liftind_CD(AR);

K1 = rho*S*Cd_min*V^3;

K2 = eff*P;

K3 = K2^2-4*W^2*V^2*Cd_min*k;

f = K1 - K2 - sqrt(K3); % Snorri (Eqn. 19-28); also see pgs. 878-879 and 883-884

end