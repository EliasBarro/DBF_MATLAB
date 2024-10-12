function [v_interp,indx] = lin_interp1(X,V,X_q,indx_st)
% v_interp,indx] = lin_interp1(X,V,X_q,indx_st)
% X are the breakpoints for V, X_q stuff is input point, indx_st is starting
% indices to make the lookup thing faster


%% find indx
indx = 0;
if nargin<4
    indx_st = 1;
end
if X(indx_st)<X_q   % search bottom-up
    for n = indx_st:length(X)         % start on the first segment
        if X(n) >= X_q                % X must be monotonically increasing
            indx = n-1;               % found the index
            break;                    % get out
        end
    end
elseif X(indx_st)>X_q   % search top-down
    for n = indx_st:-1:1              % start on the last segment
        if X(n) <= X_q                % X must be monotonically increasing
            indx = n;                 % found the index
            break;                    % get out
        end
    end
else
    indx = indx_st;
end
if indx<1
    error('X_q not found in X');
end

%% calculate the interpolation fraction
X_low   = X(indx);
X_hi    = X(indx+1);
f_x     = (X_q-X_low)/(X_hi-X_low);

% calculate the interpolated output 'v'
v_low       = V(indx);
v_hi        = V(indx+1);
v_interp    = f_x*v_hi + (1-f_x)*v_low;

end