% Author: Calina A. Copos
% Last updated: 9/14/21
clear;
%close all;
clc;

% Forward cell-cell problem
% Given: Tensions
% Unknowns: Morphology and migration speed 

% %Dimensionless parameters (asymmetric)
% f   = 0.0;%0.5;          % maximal adhesion/protrusion
% t0  = 1.5;          % characteristic tension that breaks adhesions
% tau = 0.1;          % ratio of membrane tension to adhesion/protrusion
% 
% tL  = 1.0;          % leader force
% tT  = 1.0;          % trailer force
% tcc = 0.9;          % cell-cell force
% 
% % Starting angles
% aL      = 0.429;
% aT      = 0.667;
% bL      = 0.307;
% bT      = 0.517;
% w       = 0.1;
% gamma   = 1.408;
% x       = (pi/2-gamma-w)/2;
% 
% h       = 0.2;
% 
% rL = 4.581;
% rT = 2.403;
% 
% cT = [-2,0];
%
% x2 = -0.84;%x4 - h*tan(x+gamma)
% y2 = y4 - h;

% %Dimensionless parameters (symmetric/result)
% f   = 0.1;          % maximal adhesion/protrusion
% t0  = 1.5;          % characteristic tension that breaks adhesions
% tau = 0.1;          % ratio of membrane tension to adhesion/protrusion
% 
% 
% tL  = 1.0;          % leader force
% tT  = 1.0;          % trailer force
% tcc = 0.25;%0.9;          % cell-cell force
% 
% % Starting angles
% aL      = 1.1941;
% aT      = 1.2596;
% bL      = 0.0819;
% bT      = 0.0144;
% w       = 0.1804;
% gamma   = 1.3904;
% x       = 0.3327;
% 
% h       = 0.9495;
% 
% rL = 1.510;
% rT = 1.370;
% 
% cT = [-0.4413,1.4685];

%% Dimensionless parameters (asymmetric)
f   = 0.0;          % maximal adhesion/protrusion
t0  = 1.9;          % characteristic tension that breaks adhesions
tau = 0.0;%0.1;          % ratio of membrane tension to adhesion/protrusion

tL  = 1.0;          % leader force
tT  = 1.1;          % trailer force
tcc = 0.2;%0.9;          % cell-cell force

% Starting angles
aL      = 0.712;
aT      = 0.845;
bL      = 0.050;
bT      = 0.149;
w       = 0.0;
gamma   = 1.12;
x       = (pi/2-gamma-w)/2;
h       = 0.7;

rL = 2.893;
rT = 2.152;

cT = [-2,0];

% %% Dimensionless parameters (symmetric)
% f   = 0.0;          % maximal adhesion/protrusion
% t0  = 1.9;          % characteristic tension that breaks adhesions
% tau = 0.0;%0.1;          % ratio of membrane tension to adhesion/protrusion
% 
% tL  = 1.0;          % leader force
% tT  = 1.0;          % trailer force
% tcc = 0.05;%0.9;          % cell-cell force
% 
% % Starting angles
% aL      = 0.962;
% aT      = 0.936;
% bL      = 0.477;
% bT      = 0.457;
% w       = 0.0;
% gamma   = 1.548;
% x       = (pi/2-gamma-w)/2;
% h       = 0.5;
% 
% rL = 1.581;
% rT = 1.641;
% 
% cT = [-2,0];

% Temporal discretization
dt   = 0.005;
Tend = 12.0;

% Initialization
% Set leader bubble center location
cL(1) = cT(1) + rT*cos(pi/2-bT) + rL*cos(pi/2-bL);
cL(2) = cT(2) + rT*sin(pi/2-bT) - rL*sin(pi/2-bL); 

% Get node locations    
x1 = cL(1) + rL*cos(pi/2-aL);
y1 = cL(2) + rL*sin(pi/2-aL);

x4 = cL(1) + rL*cos(pi/2+bL);
y4 = cL(2) + rL*sin(pi/2+bL);

x3 = cT(1) + rT*cos(pi/2+aT);
y3 = cT(2) + rT*sin(pi/2+aT);

x2 = x4 - h*tan(w);
%y2 = y4 - h;
%x2 = cT(1) - rT*cos(aT)*tan(bT);
%x2 = -0.84;%x4 - h*tan(x+gamma)
y2 = y4 - h;

zeta1 = 1.0;%0.05;
zeta2 = 1.0;
zeta3 = 1.0;%0.05;
zeta4 = 0.05;

% Assign result vectors from run
velx_vec    = zeros(round(Tend/dt),4);
vely_vec    = zeros(round(Tend/dt),4);
r_vec       = zeros(round(Tend/dt),3);
gamma_vec   = zeros(round(Tend/dt),1);
w_vec       = zeros(round(Tend/dt),1);

outevery = 50;

vidObj = VideoWriter('simulation_velocities_fp0dot6.mp4','MPEG-4');
open(vidObj);
for t=1:round(Tend/dt)
    
    %% Get leader cell's new morphology
    % Find radius of leader cell (nodes 1,2,4) s.t. area = 1
    clear xa ya x0 y0 a b r th da th_pre 
    da_triL = polyarea([x1,x2,x4],[y1,y2,y4]); 
    th_pre = aL+bL;
    s14 = sqrt((x1-x4)^2+(y1-y4)^2);
    for th = (th_pre-0.1):0.01:(th_pre+0.1)
        da = da_triL;
        ds = s14;
        for r = 0.1:0.01:14.99
            da = da_triL + 0.5*(th-sin(th))*r^2;
            ds = abs(s14 - 2*r*sin(th/2));
            if (da>=1 && ds<=0.01)
                %sprintf('leader: %.2f %.2f %.2f',r,th,da)
                rL = r;
                daL = da;
                break
            end
        end
        if (da>=1 && ds<=0.1)
            sprintf('leader: %.2f %.2f %.2f',r,th,da)
         	rL = r;
          	daL = da;
        	break
        end
    end

    % Find center of leader cell
    xa = 0.5*(x4-x1); ya = 0.5*(y4-y1);
    x0 = x1+xa; y0 = y1+ya;
    %x0 = x4-xa; y0 = y4-ya;
    a  = sqrt(xa^2+ya^2);
    b  = sqrt(rL^2-a^2);
    cL(1) = x0 - b*ya/a; %x0 + b*ya/a
    cL(2) = y0 + b*xa/a; %y0 - b*xa/a;
    %plot(cL(1)+rL*cos(linspace(0,2*pi,100)),cL(2)+rL*sin(linspace(0,2*pi,100)),'-.k');

    % Isolate angles alpha and beta
    aL = asin(abs(cL(1)-x1)/rL);
    bL = th-aL;

    %% Get trailer cell's new morphology
    % Find radius of trailer cell (nodes 2,3,4) s.t. area = 1
    clear xa ya x0 y0 a b r th da th_pre
    da_triT = polyarea([x3,x4,x2],[y3,y4,y2]);    
    th_pre = aT+bT;
    s34 = sqrt((x3-x4)^2+(y3-y4)^2);
    for th = (th_pre-0.1):0.01:(th_pre+0.1)
        da = da_triT;
        ds = s34;
        for r = 0.1:0.01:rL
            da = da_triT + 0.5*(th-sin(th))*r^2;
            ds = abs(s34 - 2*r*sin(th/2));
            if (da>=1 && ds<=0.01)
                %sprintf('trailer: %.2f %.2f %.2f',r,th,da)
                rT = r;
                daT = da;
                break
            end
        end
        if (da>=1 && ds<=0.01)
                sprintf('trailer: %.2f %.2f %.2f',r,th,da)
                rT = r;
                daT = da;
                break
        end
    end

    % Find center of trailer cell
    xa = 0.5*(x4-x3); ya = 0.5*(y4-y3);
    x0 = x4-xa; y0 = y4-ya;
    a  = sqrt(xa^2+ya^2);
    b  = sqrt(rT^2-a^2);
    cT(1) = x0 + b*ya/a; %x0 + b*ya/a
    cT(2) = y0 - b*xa/a; %y0 - b*xa/a;
    %plot(cT(1)+rT*cos(linspace(0,2*pi,100)),cT(2)+rT*sin(linspace(0,2*pi,100)),'-.k');
    
    % Isolate angles alpha and beta
    aT = asin(abs(x3-cT(1))/rT);
    bT = th-aT;
    
    %% Get shape at cell-cell interface
    % Find the radius at cell-cell interface
    if (tcc == 0) || ((tT/rT - tL/rL)<1e-14)
       rcc = 100;
    else
        rcc =  tcc/(tT/rT - tL/rL);
    end
    

    % Find center of circle 
    clear xa ya x0 y0 a b r th da th_pre 
    xa = 0.5*(x4-x2); ya = 0.5*(y4-y2);
    x0 = x2+xa; y0 = y2+ya;
    a  = sqrt(xa^2+ya^2);
    b  = sqrt(rcc^2-a^2);
    cC(1) = x0 - b*ya/a; %x0 + b*ya/a
    cC(2) = y0 + b*xa/a; %y0 - b*xa/a;
    %plot(cC(1)+rcc*cos(linspace(0,2*pi,100)),cC(2)+rcc*sin(linspace(0,2*pi,100)),'-.k');

    % Find angles x, gamma, and w
    s2 = (x2-x4)^2+(y2-y4)^2;
    s  = sqrt(s2);
    h  = y4-y2;
    x  = 0.5*acos(1-s2/(2*rcc*rcc)); 
    %x = asin(s/(2*rcc)); % edited 8/4
    
    if(imag(x)~=0)
        keyboard();
    end
    
    gamma = asin(h/s);
    w     = pi/2-gamma;

    %% Relabel figure and draw arcs
    th1 = linspace(pi/2-aL,pi/2+bL,100);
    %plot(cL(1)+rL*cos(th1),cL(2)+rL*sin(th1),'-','color',[0.7 0 0.5],'linewidth',2);
    th2 = linspace(pi/2+aT,pi/2-bT,100);
    %plot(cT(1)+rT*cos(th2),cT(2)+rT*sin(th2),'-','color',[0.7 0 0.5],'linewidth',2);
    th3_start = 3*pi/2+acos((cC(2)-y2)/rcc);
    th3 = linspace(th3_start,th3_start+2*x,100);
    %plot(cC(1)+rcc*cos(th3),cC(2)+rcc*sin(th3),'-','color',[0.5 0 0.5],'linewidth',2);
    %line([x1,x3],[y1,y3],'color',[0.7 0 0.5],'linewidth',2);
    
    %% Fill in shaded areas
    xLarea = [cL(1)+rL*cos(th1),cC(1)+rcc*cos(th3),linspace(x2,x1,100)];
    yLarea = [cL(2)+rL*sin(th1),cC(2)+rcc*sin(th3),y1*ones(1,100)];
    %fill(xLarea,yLarea,[0.7 0 0.5],'facealpha',0.2,'edgecolor','none');
    
    xTarea = [cT(1)+rT*cos(th2),cC(1)+rcc*cos(th3),linspace(x2,x3,100)];
    yTarea = [cT(2)+rT*sin(th2),cC(2)+rcc*sin(th3),y2*ones(1,100)];
    %fill(xTarea,yTarea,[0.7 0 0.5],'facealpha',0.2,'edgecolor','none');
    
%     scatter(x1,y1,100,'ok','fill','k')
%     scatter(x2,y2,100,'ko','fill','k')
%     scatter(x3,y3,100,'ko','fill','k')
%     scatter(x4,y4,100,'ko','fill','k')
    %title({['r_{L} = ',num2str(double(rL)),', r_{T} = ',num2str(double(rT)),', r_{cc} = ',num2str(double(rcc))],['t_{L} = ', num2str(tL), ', \alpha_{L} = ', num2str(aL),', \beta_{L} = ', num2str(bL)],['t_{T} = ', num2str(tT), ', \alpha_{T} = ', num2str(aT),', \beta_{T} = ',num2str(bT)],[', t_{cc} = ',num2str(tcc),', x = ',num2str(x),', \gamma = ', num2str(gamma),', w = ', num2str(w),', h = ', num2str(h)]});  
    
    %% Update nodes based on force balances
%     if (x2-x4)<=0
%         v2x = tcc*cos(gamma-x); % check
%         %quiver(x2,y2,tcc*cos(gamma-x),tcc*sin(gamma-x),'-m','linewidth',2);
%         
%         %v4x = tL*cos(bL) - tT*cos(bT) + tcc*sin(3*pi/2-w+x); 
%         %v4y = tcc*cos(3*pi/2-w+x) + tL*sin(bL) +  tT*sin(bT);
%         %v4x = tL*cos(bL) - tT*cos(bT) + tcc*sin(pi+gamma+x); 
%         %v4y = tcc*cos(pi+gamma+x) - tL*sin(bL) - tT*sin(bT);
%         if x<=w
%             v4x = tL*cos(bL) - tT*cos(bT) - tcc*sin(w-x); 
%             v4y = -tcc*cos(w-x) + tL*sin(bL) +  tT*sin(bT);
%             %quiver(x4,y4,-tcc*sin(w+x),-tcc*cos(w+x),'-m','linewidth',2);
%         else
%             v4x = tL*cos(bL) - tT*cos(bT) + tcc*sin(x-w); 
%             v4y = -tcc*cos(x-w) + tL*sin(bL) +  tT*sin(bT);
%             %quiver(x4,y4,tcc*sin(x-w),-tcc*cos(x-w),'-m','linewidth',2);
%         end
%     else
%        %v2x = -tcc*cos(gamma);
%        v2x = tcc*cos(gamma+x);%tcc*cos(pi-gamma-x); %check
%        %quiver(x2,y2,tcc*cos(gamma-x),tcc*sin(gamma-x),'-m','linewidth',2);
%        
%        v4x = tL*cos(bL) - tT*cos(bT) + tcc*sin(3*pi/2+x+w);
%        v4y = tcc*cos(3*pi/2+x+w) + tL*sin(bL) + tT*sin(bT);
%        %quiver(x4,y4,tcc*cos(3*pi/2+w+x),tcc*sin(3*pi/2+w+x),'-m','linewidth',2);
%     end
    
    if (x2-x4)<=0
        v2x = tcc*cos(gamma-x);
        v4x = tL*cos(bL) - tT*cos(bT) - tcc*sin(w-x);
        v4y = -tcc*cos(w-x) + tL*sin(bL) +  tT*sin(bT);
        
        %quiver(x2,y2,tcc*cos(gamma-x),tcc*sin(gamma-x),'-k','linewidth',1);
        %quiver(x4,y4,-tcc*sin(w-x),-tcc*cos(w-x),'-k','linewidth',1);
    else
        v2x = -tcc*cos(gamma+x);
        v4x = tL*cos(bL) - tT*cos(bT) + tcc*sin(w+x);
        v4y = -tcc*cos(w+x) + tL*sin(bL) +  tT*sin(bT);
        
        %quiver(x2,y2,-tcc*cos(gamma+x),tcc*sin(gamma+x),'-k','linewidth',1);
        %quiver(x4,y4,tcc*sin(w+x),-tcc*cos(w+x),'-k','linewidth',1);
    end     

    %v1x = (1+f)*(1-tL/(t0+tau)) - tL*cos(aL); % check
    v1x = 0.6 - tL*cos(aL); % check
    v1x = v1x/zeta1;
    %v3x = tT*cos(aT) - (1-f)*(1-tT/(t0-tau)); % check
    v3x = tT*cos(aT) - 0.0; % check
    v3x = v3x/zeta3;
    
    v2x = v2x/zeta2;
    v4x = v4x/zeta4;
    v4y = v4y/zeta4;
    
    %quiver(x4,y4,tL*cos(bL),tL*sin(bL),'-m','linewidth',2);
    %quiver(x4,y4,-tT*cos(bT),tT*sin(bT),'-m','linewidth',2);
    
%     % forces
%     quiver(x1,y1,0.3,0,'-k','linewidth',1);
%     quiver(x1,y1,-tL*cos(aL),tL*sin(aL),'-k','linewidth',1);
%     quiver(x3,y3,tT*cos(aT),tT*sin(aT),'-k','linewidth',1);
%     quiver(x4,y4,tL*cos(bL),tL*sin(bL),'-k','linewidth',1);
%     quiver(x4,y4,-tT*cos(bT),tT*sin(bT),'-k','linewidth',1)
%    
    
    if mod(t,outevery) == 0
        %% Plot nodes and `vertex model schematic'
        figure(1);
        scatter(x1,y1,100,'ok','fill','k')
        hold on;
        scatter(x2,y2,100,'ko','fill','k')
        scatter(x3,y3,100,'ko','fill','k')
        scatter(x4,y4,100,'ko','fill','k')
        set(gcf,'color','w');
        set(gca,'fontsize',20,'fontname','Times New Roman'); box on;
        %line([x1,x2],[y1,y2],'color','b','linewidth',1); line([x2,x3],[y2,y3],'color','b','linewidth',1); line([x3,x4],[y3,y4],'color','b','linewidth',1); line([x4,x1],[y4,y1],'color','b','linewidth',1); line([x4,x2],[y4,y2],'color','b','linewidth',1);
        xlim([-4 4]); ylim([1 3]);
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        pbaspect([8 2 1]);
        
        plot(cL(1)+rL*cos(th1),cL(2)+rL*sin(th1),'-','color',[0.7 0 0.5],'linewidth',2);
        plot(cT(1)+rT*cos(th2),cT(2)+rT*sin(th2),'-','color',[0.7 0 0.5],'linewidth',2);
        plot(cC(1)+rcc*cos(th3),cC(2)+rcc*sin(th3),'-','color',[0.5 0 0.5],'linewidth',2);
        line([x1,x3],[y1,y3],'color',[0.7 0 0.5],'linewidth',2);

        fill(xLarea,yLarea,[0.7 0 0.5],'facealpha',0.4,'edgecolor','none');
        fill(xTarea,yTarea,[0.7 0 0.5],'facealpha',0.2,'edgecolor','none');

        scatter(x1,y1,100,'ok','fill','k')
        scatter(x2,y2,100,'ko','fill','k')
        scatter(x3,y3,100,'ko','fill','k')
        scatter(x4,y4,100,'ko','fill','k')

        % velocities
        quiver(x1,y1,v1x,0,'-','color',[0.75 0.75 0.75],'linewidth',3,'maxheadsize',10);
        quiver(x2,y2,v2x,0,'-','color',[0.75 0.75 0.75],'linewidth',3,'maxheadsize',10);
        quiver(x3,y3,v3x,0,'-','color',[0.75 0.75 0.75],'linewidth',3,'maxheadsize',10);
        %quiver(x4,y4,v4x,v4y,'-','color',[54/256 69/256 79/256],'linewidth',3,'maxheadsize',10);
        plot(1*ones(100,1),linspace(0,4,100),'-','color',[0.5 0.5 0.5],'linewidth',0.1);
        currFrame = getframe(gcf);
        hold off;
        writeVideo(vidObj,currFrame);
    end
    
    
    x1  = x1 + dt*v1x;
    x2  = x2 + dt*v2x;
    x3  = x3 + dt*v3x;
    x4  = x4 + dt*v4x;
    y4  = y4 + dt*v4y;
    
    velx_vec(t,:) = [v1x,v2x,v3x,v4x];
    vely_vec(t,:) = [0,0,0,v4y];
    r_vec(t,:) = [rL,rT,rcc];
    gamma_vec(t) = gamma;
    w_vec(t) = w;
end
close(vidObj);