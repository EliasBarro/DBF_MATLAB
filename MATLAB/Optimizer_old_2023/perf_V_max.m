% Michael Qiu - 2022

function [V_max, counter] = perf_V_max(ALT,S,Cd_min,AR,W,P,eff)

V0 = 0;
V1 = 50;
f0 = -1000*eff*P;

counter = 0;

while abs(V1 - V0) > 0.0001 && counter < 1000
    
    Vmid = 0.5*(V0+V1);
    fmid = perf_f_of_V(ALT,S,Cd_min,AR,W,P,eff,Vmid);
    
    if f0*fmid < 0
        V1 = Vmid;
    else
        V0 = Vmid;
        f0=fmid;
    end
    
    counter = counter + 1;
    
end

V_max = V1;

end