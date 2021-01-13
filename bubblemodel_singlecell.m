clear;
close all;
clc;

% Model parameters
kappa1 = 1;
kappa2 = 1;
kappa3 = 1;
kappa4 = 1;
T = 1000;% cortical tension (pN/um)
V = 1000;% cell volume (um^3)

dt = 0.01;
Tmax = 1;
Nt = Tmax/dt;

r = 10;
th = 0:0.01:2*pi;
X0 = r*cos(th) + 5;
Y0 = r*sin(th) + 5;
X = X0; Y = Y0;

for i = 1:Nt
    i
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Setup Poisson's equation with Dirichlet BC
    % to solve the height of the droplet and its pressure
    model = createpde();
    
    % create geometry
    %geometryFromEdges(model,@circleg);
    % TBD (done? 7/1): figure out how to put boundary points only and create
    % geometry from those boundary points
    pgon = polyshape(X,Y);
    tr = triangulation(pgon);
    tnodes = tr.Points';
    telements = tr.ConnectivityList';
    geometryFromMesh(model,tnodes,telements);
    
    %figure
    %pdegplot(model,'EdgeLabels','on');
    %axis equal

    % setup BC and PDE coefficients
    applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0);
    specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',1);

    % create mesh 
    hmax = 0.5;
    mesh = generateMesh(model,'Hmax',hmax,'GeometricOrder','linear');
    %figure
    %pdemesh(model);
    %axis equal

    % solve PDE
    results = solvepde(model);
    u = results.NodalSolution;
    figure(1); 
    pdeplot(model,'XYData',u);
    hold on;
    title('H = h(x(t),y(t))P(t)/T');
    xlabel('x');
    ylabel('y');
    hold off;

    % get current area
    [a,aelements] = area(mesh);
    
    % distribute area per node
    anodes = zeros(length(mesh.Nodes),1);
    for j = 1:length(mesh.Elements)
        anodes(mesh.Elements(j)) = anodes(mesh.Elements(j)) + 1/3*aelements(j);
    end
    
    % get initial area
    if i==1
        a0 = a;
    end

    % find pressure in droplet
    pressure = T*sum(u.*anodes)/V;

    % find contact angle around the boundary 
    % TBD: boundary only!!!!
    
    bdd = convhull(mesh.Nodes(1,:),mesh.Nodes(2,:));
    X = mesh.Nodes(1,bdd);
    Y = mesh.Nodes(2,bdd);
    querypoints = [X;Y];
    [gradx,grady] = evaluateGradient(results,querypoints);
    phi = atan2(grady,gradx);
    
    % find center of mesh
    xc = mean(mesh.Nodes(1,:));
    yc = mean(mesh.Nodes(2,:));
    costh = (X-xc)./sqrt((X-xc).^2+(Y-yc).^2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. Update velocity of droplet in the normal direction
    v = kappa1*(a0-a) + kappa2*costh - kappa3*T*cos(phi)' + kappa4*pressure;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. Update the boundary using a level-set function
    % L_t + v |grad L| = 0.
    % TBD: How to do this?
    [XX,YY] = meshgrid(X,Y);
    L = sqrt(XX.^2 + YY.^2);
    L(L>100) = 100;
    L(L<99) = 99;
    L = L - 99.5;
    LS = levelset2D(L, 20);
    LS = reinitialize(LS);
    figure(2);
    plot(LS, 'contour', 'levelset');
    %LS = propagate(LS,1,'speed_normal', v);
    
end




