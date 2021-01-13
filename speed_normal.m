function [dphi_dt,dt] = speed_normal(ls, varargin)
    
    % Assume that the speed is passed as argument
    F = varargin{1};
    if isscalar(F) % expand if given as a scalar
        F = F*ones(size(ls));
    end

    % Determine stable timestep for explicit time integration
    dt = 0.9/max(abs(F(:)));

    % Compute upwind differences using Godunov's method
    [Dx2,Dy2] = godunov(ls,F);

    % Compute and return level set PDE (dphi_dt)
    dphi_dt = -F .* sqrt(Dx2 + Dy2);
end