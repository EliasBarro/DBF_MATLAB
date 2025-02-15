function W = reqpower_mod(TW,AUW,WL,CLmax,e_eff,p_eff,ALT)

g = 9.81;

VLOF = 1.1*sqrt((2*WL)/(dens(ALT)*CLmax));

W = (TW*(AUW*g)).*(VLOF/sqrt(2))./(e_eff*p_eff); %Modified from Snorri pg. 59

end