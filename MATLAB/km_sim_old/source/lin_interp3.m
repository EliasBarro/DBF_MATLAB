function [v_interp,indx,indy,indz] = lin_interp3(X,Y,Z,V, X_q,Y_q,Z_q, indx_st,indy_st,indz_st)

% x y and z are breakpoints for V, _q stuff is input point, _st is starting
% indices to make the lookup thing faster

%% find indices
% find indx
indx = 0;
if nargin<8
    indx_st = 1;
end
if X(indx_st)<X_q   % search bottom-up
    for n = indx_st:length(X)         % start on the first segment
        if X(n) >= X_q                  % X must be monotonically increasing
            indx = n-1;               % found the index
            break;                    % get out
        end
    end
elseif X(indx_st)>X_q   % search top-down
    for n = indx_st:-1:1              % start on the last segment
        if X(n) <= X_q                % X must be monotonically increasing
            indx = n;               % found the index
            break;                    % get out
        end
    end
else
    indx = indx_st;
end
if indx<1
    error('X_q not found in X');
end

% find indy
indy = 0;
if nargin<9
    indy_st = 1;
end
if Y(indy_st)<Y_q   % search bottom-up
    for n = indy_st:length(Y)         % start on the first segment
        if Y(n) >= Y_q                % Y must be monotonically increasing
            indy = n-1;               % found the index
            break;                    % get out
        end
    end
elseif Y(indy_st)>Y_q   % search top-down
    for n = indy_st:-1:1              % start on the last segment
        if Y(n) <= Y_q                 % Y must be monotonically increasing
            indy = n;               % found the index
            break;                    % get out
        end
    end
else
    indy = indy_st;
end
if indy<1
    error('Y_q not found in Y (plane is flying too fast or too slow)');
end

% find indz
indz = 0;
if nargin<10
    indz_st = 1;
end
if Z(indz_st)<Z_q   % search bottom-up
    for n = indz_st:length(Z)         % start on the first segment
        if Z(n) >= Z_q                % Z must be monotonically increasing
            indz = n-1;               % found the index
            break;                    % get out
        end
    end
elseif Z(indz_st)>Z_q   % search top-down
    for n = indz_st:-1:1              % start on the last segment
        if Z(n) <= Z_q                 % Z must be monotonically increasing
            indz = n;               % found the index
            break;                    % get out
        end
    end
else
    indz = indz_st;
end
if indz<1
    error('Z_q not found in Z');
end

%% calculate interpolation fractions
% f_x
X_lo    = X(indx);
X_hi    = X(indx+1);
f_x     = (X_q-X_lo)/(X_hi-X_lo);

% f_y
Y_lo    = Y(indy);
Y_hi    = Y(indy+1);
f_y     = (Y_q-Y_lo)/(Y_hi-Y_lo);

% f_z
Z_lo    = Z(indz);
Z_hi    = Z(indz+1);
f_z     = (Z_q-Z_lo)/(Z_hi-Z_lo);

%% define bounding points and calculate interpolated value
x0y0z0 = V(indx+0,indy+0,indz+0);
x0y0z1 = V(indx+0,indy+0,indz+1);
x0y1z0 = V(indx+0,indy+1,indz+0);
x0y1z1 = V(indx+0,indy+1,indz+1);
x1y0z0 = V(indx+1,indy+0,indz+0);
x1y0z1 = V(indx+1,indy+0,indz+1);
x1y1z0 = V(indx+1,indy+1,indz+0);
x1y1z1 = V(indx+1,indy+1,indz+1);

v_interp = (1-f_x)*(1-f_y)*(1-f_z)*x0y0z0 + ...
           (1-f_x)*(1-f_y)*(  f_z)*x0y0z1 + ...
           (1-f_x)*(  f_y)*(1-f_z)*x0y1z0 + ...
           (1-f_x)*(  f_y)*(  f_z)*x0y1z1 + ...
           (  f_x)*(1-f_y)*(1-f_z)*x1y0z0 + ...
           (  f_x)*(1-f_y)*(  f_z)*x1y0z1 + ...
           (  f_x)*(  f_y)*(1-f_z)*x1y1z0 + ...
           (  f_x)*(  f_y)*(  f_z)*x1y1z1;

end