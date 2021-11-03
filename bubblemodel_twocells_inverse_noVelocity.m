%% Author: Calina A. Copos
% Last updated: 7/8/21
clear;
close all;
clc;

% Dimensionless parameters
f   = 0.2;          % maximal adhesion/protrusion
t0  = 1.5;          % characteristic tension that breaks adhesions
tau = 0.1;          % ratio of membrane tension to adhesion/protrusion

%% Numerical solver
t1  = 1.0;          % leader force
t2  = 1.0;          % trailer force
tcc = 0.7;          % cell-cell force

h   = 0.6;
w   = 0.2;

% Step 1. Force balance at the top of cell-cell region.
% (1.a') t1 cos(phi1) - t2 cos(phi2) = 0
% (1.b') t1 sin(phi1) + t2 sin(phi2) - tcc = 0
%func1 = @(phi2) t1*sin( acos((t2/t1)*cos(phi2)) ) + t1*sin(phi2) - tcc;
%phi2 = fzero(func1,[-pi/2,pi/2]);
%phi1 = acos( (t2/t1)*cos(phi2) );

% (1.a'') x = cos(phi1), y = cos(phi2), x = t2/t1 * y
% (1.b'') sqrt(1-(t2/t1)^2*y^2) + (t2/t1)*sqrt(1-y^2) - tcc = 0
%func1 = @(y) sqrt(1-((t2/t1)^2)*y^2) + (t2/t1)*sqrt(1-y^2) - tcc;
%y = fzero(func1,0.8);
%phi2 = acos(y);
%phi1 = acos((t2/t1)*y);

% alternative because fzero is not working
syms y;
func1 = sqrt(1-((t2/t1)^2)*y^2) + (t2/t1)*sqrt(1-y^2) - tcc;
ys = vpasolve(func1,y,[0 Inf]);
phi2 = acos(ys);
phi1 = acos((t2/t1)*ys);

b1 = phi1 + w;
b2 = phi2 - w;

% Step 2. Trigonometry & Step 3. Area conservation
% (2.a) h = h1 = r1*(cos(b1)-cos(a1))
% (2.b) h = h2 = r2*(cos(b2)-cos(a2))
% (2.c) area of leader - 1 = da1 = 0
% (2.d) area of trailer - 1 = da2 = 0 

%da1 = @(r1) r1^2*( acos(cos(b1)-h/r1) + b1 ) - ...
%            r1^2*sin(acos(cos(b1)-h/r1))*(cos(b1)-h/r1) - ...
%            r1^2*(cos(b1)-h/r1)^2*tan(b1) + ...
%            r1^2*( cos(b1)-(cos(b1)-h/r1)*tan(b1) )*( sin(b1)-(cos(b1)-h/r1)*tan(b1)*tan(b1) ) - ...
%            2;
%try 
%    r1 = fzero(da1,3);
%catch 
%    continue
%end 

% alternative because fzero is not working
syms rs1;
%dA1 = rs1^2*( acos(cos(b1)-h/rs1) + b1 ) - rs1^2*sin(acos(cos(b1)-h/rs1))*(cos(b1)-h/rs1) - rs1^2*(cos(b1)-h/rs1)^2*tan(b1) + rs1^2*( cos(b1)-(cos(b1)-h/rs1)*tan(b1) )*( sin(b1)-(cos(b1)-h/rs1)*tan(b1)*tan(b1) ) - 2;
dA1 = rs1^2*( acos(cos(b1)-h/rs1) + b1 ) - rs1^2*sin(acos(cos(b1)-h/rs1))*(cos(b1)-h/rs1) - rs1^2*(cos(b1)-h/rs1)^2*tan(b1) - 2;
r1 = vpasolve(dA1,rs1,[0 Inf]);
a1  = acos(cos(b1) - h/r1);

%da2 = @(r2) r2^2*( acos(cos(b2)-h/r2) + b2 ) - ...
%            r2^2*sin(acos(cos(b2)-h/r2))*(cos(b2)-h/r2) - ...
%            r2^2*(cos(b2)-h/r2)^2*tan(b2) + ...
%            r2^2*( cos(b2)-(cos(b2)-h/r2)*tan(b2) )*( sin(b2)-(cos(b2)-h/r2)*tan(b2)*tan(b2) ) - ...
%            2;
%try 
%    r2 = fzero(da2,1);
%catch 
%    continue
%end

% alternative because fzero is not working
syms rs2;
%dA2 = rs2^2*( acos(cos(b2)-h/rs2) + b2 ) - rs2^2*sin(acos(cos(b2)-h/rs2))*(cos(b2)-h/rs2) -rs2^2*(cos(b2)-h/rs2)^2*tan(b2) + rs2^2*( cos(b2)-(cos(b2)-h/rs2)*tan(b2) )*( sin(b2)-(cos(b2)-h/rs2)*tan(b2)*tan(b2) ) - 2;
dA2 = rs2^2*( acos(cos(b2)-h/rs2) + b2 ) - rs2^2*sin(acos(cos(b2)-h/rs2))*(cos(b2)-h/rs2) -rs2^2*(cos(b2)-h/rs2)^2*tan(b2) - 2;
r2 = vpasolve(dA2,rs2,[0 Inf]);
a2 = acos(cos(b2) - h/r2);

% Step 3. Find the radius of curvature for the cell-cell junction:
if (tcc == 0) || (abs(t2/r2 - t1/r1)<1e-14)
   rcc = 100;
else
    rcc =  tcc/abs(t2/r2 - t1/r1);
end

% Step 4. From trigonometry at the cell-cell interface   
%func2 = @(x) h^2/(cos(abs(x)+abs(w))^2) - 2*rcc^2*(1-cos(2*abs(x)));    
%try 
%    x = fzero(func2,0);
%catch 
%    continue
%end

% alternative because fzero is not working
syms ex;
func2 = h^2/(cos(ex+abs(w))^2) - 2*rcc^2*(1-cos(2*ex));
x = vpasolve(func2,ex,[0,pi/4]);
gamma = pi/2-abs(w)-2*x;

% Step 5. Force balance at the ventral points.
zeta1   = 0.75;          % effective drag on leader cell
zeta2   = 0.75;          % effective drag on trailer cell
zetacc  = 0.75;          % effective drag at cell-cell zone
%f1 = @(v) (1+f)*(1-t1/(t0+tau)) - zeta1*v - t1*cos(a1)
%f2 = @(v) (1-f)*(1-t2/(t0-tau)) + zeta2*v - t2*cos(a2)
%fcc= @(v) tcc*cos(gamma) - zetacc*v

%v1 = fzero(f1,0);
%v2 = fzero(f2,0);
%vcc= fzero(fcc,0);

v1 = ((1+f)*(1-t1/(t0+tau)) - t1*cos(a1));
v2 = (t2*cos(a2) - (1-f)*(1-t2/(t0-tau)));
v3 = tcc*cos(gamma)/zetacc;

%% Plotting
wsgn = w/abs(w);
% Trailer cell
%center2 = [-2+r2*cos(a2),0-r2*sin(a2)];
center2 = [-2,0];
start_angle2 = pi/2+a2;
end_angle2 = pi/2-b2;
N = 100;
th2 = linspace(start_angle2,end_angle2,N);
x2 = center2(1) + r2*cos(th2);
y2 = center2(2) + r2*sin(th2);
plot(x2,y2,'-','color',[0.7 0 0.5],'linewidth',2);
hold on

% Leader cell
center1 = [center2(1)+r2*cos(pi/2-b2)+r1*cos(pi/2-b1);center2(2)-(r1*sin(pi/2-b1)-r2*sin(pi/2-b2))];
start_angle1 = pi/2-a1;
end_angle1 = pi/2+b1;
th1 = linspace(start_angle1,end_angle1,N);
x1 = center1(1) + r1*cos(th1);
y1 = center1(2) + r1*sin(th1);
plot(x1,y1,'-','color',[0.7 0 0.5],'linewidth',2);
hold on;

if (w==0 || wsgn<0)
    start_angle_cc = 3*pi/2+gamma+2*x;
    end_angle_cc = start_angle_cc-2*x;
elseif wsgn>0
    start_angle_cc = 3*pi/2-gamma-2*x;
    end_angle_cc = start_angle_cc+2*x;
end

start_angle_cc = 3*pi/2+gamma+2*abs(x);
end_angle_cc = start_angle_cc-2*abs(x);

th3 = linspace(start_angle_cc,end_angle_cc,N);
center3 = [x1(end)-h*cot(gamma+abs(x))-rcc*cos(pi/2-gamma);y1(end)+rcc*sin(abs(w))];
x3 = center3(1) + rcc*cos(th3);
y3 = center3(2) + rcc*sin(th3);
plot(x3,y3,'-','color',[0.7 0 0.5],'linewidth',2);

x1area = [x1,x3];
y1area = [y1,y3];
fill(x1area,y1area,[0.7 0 0.5],'facealpha',0.4,'edgecolor','none');

x2area = [x2,x3];
y2area = [y2,y3];
fill(x2area,y2area,[0.7 0 0.5],'facealpha',0.2,'edgecolor','none');

plot(linspace(-4,6,100),y2(1)*ones(100,1),'-','color',[0.75 0.75 0.75],'linewidth',3);
scatter(x1(1),y1(1),100,'ko','fill');
scatter(x2(1),y2(1),100,'ko','fill');
scatter(x1(end),y1(end),100,'ko','fill');
set(gca,'fontsize',20,'fontname','Times New Roman'); box on;
set(gcf,'color','w');
%title({['\zeta_{cc} = ',num2str(zetacc), ', r_{1} = ',num2str(r1),', r_{2} = ',num2str(r2),', r_{cc} = ',num2str(rcc)],['t_{1} = ', num2str(t1), ', \alpha_{1} = ', num2str(a1*180/pi),', \beta_{1} = ', num2str(real(b1)*180/pi)],['t_{2} = ', num2str(t2), ', \alpha_{2} = ', num2str(a2*180/pi),', \beta_{2} = ',num2str(real(b2)*180/pi)],['\gamma = ', num2str(real(gamma)*180/pi),', w = ', num2str(real(w)*180/pi),', v1 = ', num2str(v1),', v2 = ', num2str(v1),', v3 = ', num2str(v3),', h = ', num2str(h),', t_{cc} = ',num2str(tcc)]});
title({['t_{L} = ', num2str(t1), ', t_{T} = ', num2str(t2), ', t_{cc} = ', num2str(tcc)],['r_{L} = ',num2str(double(r1)),', r_{T} = ',num2str(double(r2)),', r_{cc} = ',num2str(double(rcc))],['\alpha_{L} = ', num2str(double(a1)),', \beta_{L} = ', num2str(double(b1)),', \alpha_{T} = ', num2str(double(a2)),', \beta_{T} = ',num2str(double(b2))],['\gamma = ', num2str(double(gamma)),', w = ', num2str(double(w)),', h = ', num2str(h)]});
xlim([-4 4]); ylim([-0.5 2.5]);
set(gca,'xtick',[])
set(gca,'ytick',[])
pbaspect([8 2 1]);

% Compute area of trailer cell (updated 6/16)
xT = [x2,x3,flip(x2)];
yT = [y2,y3,y2(1)*ones(1,length(y2))];
aT = polyarea(xT,yT);

% Compute area of leader cell (updated 6/16)
xL = [x1,x3,flip(x1)];
yL = [y1,y3,y1(1)*ones(1,length(y1))];
aL = polyarea(xL,yL);

sprintf('Area of leader cell %.4f and trailer cell %.4f',aL,aT)
