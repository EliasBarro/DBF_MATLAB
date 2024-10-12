function T_s = static_thrust(v_c, T_c, v_lof, T_lof)
    fun = @(Ts,Tv,v) (Ts - 2*T_c)./(v_c.^2)*v.^2 + (3*T_c - 2*Ts)./v_c*v + Ts - Tv;
    T_s = fzero(@(x) fun(x,T_lof,v_lof),40);
    fcn = @(v) (T_s - 2*T_c)./(v_c.^2)*v.^2 + (3*T_c - 2*T_s)./v_c*v + T_s;
end