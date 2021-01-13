%% Level set operator defining a PDE for advection 
% Implements the level set PDE for advection used in the example:
% <advect3D.html advect3D.m>

function [dphi_dt,dt] = advect3D_operator(ls, varargin)

% Assume that the components of the external vector field is passed as two
% arguments
Fx = varargin{1};
Fy = varargin{2};
Fz = varargin{3};

% Extract components only within narrowband
Fx = Fx(narrowband(ls));
Fy = Fy(narrowband(ls));
Fz = Fz(narrowband(ls));

% Determine stable timestep for explicit time integration
normgrad = sqrt(Fx.^2 + Fy.^2 + Fz.^2);
dt = 0.9 / max(normgrad(:));

% Compute differences for upwinding
[Dx_m,Dx_p,Dy_m,Dy_p,Dz_m,Dz_p] = diff_upwind(ls);

% Determine correct difference according to upwind direction (x component)
mask = Fx<0;
Dx = zeros(size(narrowband(ls)));
Dx(mask) = Dx_p(mask);
Dx(~mask) = Dx_m(~mask);

% Determine correct difference according to upwind direction (y component)
mask = Fy<0;
Dy = zeros(size(narrowband(ls)));
Dy(mask) = Dy_p(mask);
Dy(~mask) = Dy_m(~mask);

% Determine correct difference according to upwind direction (z component)
mask = Fz<0;
Dz = zeros(size(narrowband(ls)));
Dz(mask) = Dz_p(mask);
Dz(~mask) = Dz_m(~mask);

% Compute and return level set PDE (dphi_dt)
dphi_dt = -(Dx.*Fx + Dy.*Fy + Dz.*Fz);
