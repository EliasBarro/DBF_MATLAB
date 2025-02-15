function W = reqpower(TW,AUW,V,e_eff,p_eff)

g = 9.81;

W = (TW*(AUW*g))*V./(e_eff*p_eff); %Modified from Snorri pg. 59

end