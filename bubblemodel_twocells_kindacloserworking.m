% Author: Calina A. Copos
% Last updated: 1/5/21
clear;
close all;
clc;

% Known parameters: f, xi, t1, t2, tcc, t0, tau

% Dimensionless parameters
f   = 0.1;          % maximal adhesion/protrusion
xi  = 0.1;          % effective drag
t1  = 1.0;          % leader force
t2  = 1.1;          % trailer force
tcc = 0.3;          % cell-cell force
t0  = 1.3;          % characteristic tension that breaks adhesions
tau = 0.1;          % ratio of membrane tension to adhesion/protrusion

epsilon = 0.001;
v        = f + t1*(tau/t0^2-f/t0);       % velocity (linear approx. v = f + t*(tau/t0^2 - f/t0))
dh       = 1.;
dv       = 1.;
dg       = 1.;
count    = 0;

%% Numerical solver

tol = 0.02; 
w = -0.5;
v = f + t1*(tau/t0^2-f/t0);
%v = 0.1;
while(dv>=tol && count<=10000)
%while(count<=1)
    dh = 1.;
    %v = 0.001;
    v = v + epsilon;
    w = -0.5;
    %w = w + epsilon;
    %v = epsilon;
    
    while(dh>tol && count<=10000)
        %v = v + epsilon;
        w = w + epsilon;
        wsgn = w/abs(w);
        
        % Step 1. Guess v. Force balance at front and rear endpoints
        % (1.a) (1+f) (1 - t1/(t0+tau)) - v = t1 cos(a1)
        % (1.b) (1-f) (1 - t2/(t0-tau)) + v = t2 cos(a2)

        a1 = acos( (1+f)*(1-t1/(t0+tau))/t1 - v/t1 );
        a2 = acos( (1-f)*(1-t2/(t0-tau))/t2 + v/t2 );
        % varies greatly with t0 (critical adhesion breaking force)

        % Step 2. Guess w. Force balance at the top of cell-cell region.
        % (2.a) t1 cos(b1) - t2 cos(b2) +/- tcc sin(|w|) = 0
        % (2.b) t1 sin(b1) + t2 sin(b2) -   tcc cos(|w|) = 0
        % play around with tcc, xi
        %func1 = @(b2) t1*sin( acos(t2/t1*cos(b2) + tcc/t1*sin(w)) ) + t2*sin(b2) - tcc*cos(w);
        %b2 = fzero(func1,[-pi/2,pi/2]);
        %b1 = acos( t2/t1*cos(b2) + tcc/t1*sin(w) );
        
        func1 = @(b2) t1*sin( acos(t2/t1*cos(b2) - wsgn*tcc/t1*sin(abs(w))) ) + t2*sin(b2) - tcc*cos(abs(w));
        %b2 = fzero(func1,[-pi/2,pi/2]);
        try 
            b2 = fzero(func1,0);
        catch 
            continue
        end
        b1 = acos( t2/t1*cos(b2) - wsgn*tcc/t1*sin(abs(w)) );
        
%         if (b1<0 || b2<0 || imag(b1)>1e-10 || imag(b2)>1e-10)
%             wsgn = -1;
%             func1 = @(b2) t1*sin( acos(t2/t1*cos(b2) - tcc/t1*sin(w)) ) + t2*sin(b2) - tcc*cos(w);
%             b2 = fzero(func1,[-pi/2,pi/2]);
%             b1 = acos( t2/t1*cos(b2) - tcc/t1*sin(w) );
%             sprintf('here')
%         end
        
        % Zero out the imaginary components
        if imag(b1)<=1e-10
            b1 = real(b1);
        end
        if imag(b2)<=1e-10
            b2 = real(b2);
        end

        % For vertical tcc (aka w approx 0)
        % (2.0.a) t1 cos(b1) = t2 cos(b2) -> b1 = acos( t2/t1 cos(b2) )
        % (2.0.b) t1 sin(b1) + t2 sin(b2) = tcc
    %     func = @(b2) t1*sin(acos(t2*cos(b2)/t1)) + t2*sin(b2) - tcc;
    %     b2 = fzero(func,-1.5);
    %     b1 = acos(t2/t1*cos(b2));

        % Step 3.
        % (3.a) r1 = (1/2 - 1/4 sin(2*a1) + (1/2 - 1/4 sin (2*b1))^(-1/2)
        % (3.b) r2 = (1/2 - 1/4 sin(2*a2) + (1/2 - 1/4 sin (2*b2))^(-1/2)
        
        %r1 = ( (0.5 - 0.25*sin(2*a1)) + (0.5 - 0.25*sin(2*b1)) )^(-1/2);
        %r2 = ( (0.5 - 0.25*sin(2*a2)) + (0.5 - 0.25*sin(2*b2)) )^(-1/2);
        
        % (attempt 6/8/21 did not work also not correct)
        %r1 = ( (a1 - 0.25*sin(2*a1)) + (b1 - 0.25*sin(2*b1)) )^(-1/2);
        %r2 = ( (a2 - 0.25*sin(2*a2)) + (b2 - 0.25*sin(2*b2)) )^(-1/2);
        
        % WRONG (re-visited)
        % (3.a) r1 = ( (a1+b1) - 0.5*sin(a1)*cos(a1) - 0.5*cos(a1)*cos(a1)*tan(b1) )^(-1/2)
        % (3.b) r2 = ( (a2+b2) - 0.5*sin(a2)*cos(a2) - 0.5*cos(a2)*cos(a2)*tan(b2) )^(-1/2)
        
        % (works 6/9/21)
        r1 = ( (a1+b1) - 0.5*sin(a1)*cos(a1) - 0.5*cos(a1)*cos(a1)*tan(b1) )^(-1/2);
        r2 = ( (a2+b2) - 0.5*sin(a2)*cos(a2) - 0.5*cos(a2)*cos(a2)*tan(b2) )^(-1/2);
        
        % Step 4.
        % (4.a) h = h1 = r1*(cos(b1)-cos(a1))
        % (4.b) h = h2 = r2*(cos(b2)-cos(a2))

        h1 = r1*(cos(b1)-cos(a1));
        h2 = r2*(cos(b2)-cos(a2));
        dh = abs(h2-h1);
        
        if(real(a1*180/pi)>=170)
            keyboard();
        end

        sprintf('Angles (degrees): a1=%.4f a2=%.4f b1=%.4f b2=%.4f w=%.4f',a1*180/pi,a2*180/pi,b1*180/pi,b2*180/pi,w*180/pi)
        sprintf('Cell-cell region: h1=%.4f h2=%.4f dh=%.4f',h1,h2,dh)
    end
    
    h = h2;

    % Step 5. Intracellular pressure
    p1 = t1/r1;
    p2 = t2/r2;

    % Step 6. Find the radius of curvature for the cell-cell junction:
    rcc = tcc/abs(t2/r2 - t1/r1);

    % Step 7.
    func2 = @(x) h^2/(cos(x+abs(w))^2) - 2*rcc^2*(1-cos(2*x));
    %func2 = @(x) h^2/(cos(x-w)^2) - 2*rcc^2*(1-cos(2*x));
    
    try 
        x = fzero(func2,0);
    catch 
        continue
    end
  
    %x = fzero(func2,[0,pi/4]);
    %gamma = pi/2-w-2*x;
    gamma = pi/2-abs(w)-2*x;
    
%     if(wsgn<0)
%         func2 = @(x) h^2/(cos(x-w)^2) - 2*rcc^2*(1-cos(2*x));
%         x = fzero(func2,[0,pi/4]);
%         gamma = pi/2-2*x+w;
%     end

    % Step 8.
    % tcc cos(gamma) = xi*v at the bottom of the cell-cell boundary
    vv = tcc*cos(gamma)/xi;
    dv = abs(vv-v);
    
    sprintf('Speed: v=%.4f vv=%.4f',v,vv)
    sprintf('Angle: gamma=%.4f',gamma)
    count = count + 1;
end

%% Plotting
% Trailer cell
%center2 = [-2+r2*cos(a2),0-r2*sin(a2)];
center2 = [-2,0];
start_angle2 = pi/2+a2;
end_angle2 = pi/2-b2;
N = 100;
th2 = linspace(start_angle2,end_angle2,N);
x2 = center2(1) + r2*cos(th2);
y2 = center2(2) + r2*sin(th2);
plot(x2,y2,'-b','linewidth',2);
hold on
%plot(x2(end)*ones(N,1),linspace(y2(1),y2(end),N),'-k','linewidth',2);

% Leader cell
center1 = [center2(1)+r2*cos(pi/2-b2)+r1*cos(pi/2-b1);center2(2)-(r1*sin(pi/2-b1)-r2*sin(pi/2-b2))];
start_angle1 = pi/2-a1;
end_angle1 = pi/2+b1;
th1 = linspace(start_angle1,end_angle1,N);
x1 = center1(1) + r1*cos(th1);
y1 = center1(2) + r1*sin(th1);
plot(x1,y1,'-b','linewidth',2);
hold on;
%plot(x1(end)*ones(N,1),linspace(y1(1),y1(end),N),'-k','linewidth',2);

%start_angle_cc = 3*pi/2+gamma;
if wsgn<0
    start_angle_cc = 3*pi/2+gamma+2*x;
    end_angle_cc = start_angle_cc-2*x;
elseif wsgn>0
    start_angle_cc = 3*pi/2-gamma-2*x;
    end_angle_cc = start_angle_cc+2*x;
end

start_angle_cc = 3*pi/2+gamma+2*x;
end_angle_cc = start_angle_cc-2*x;

th3 = linspace(start_angle_cc,end_angle_cc,N);
%center3 = [x1(end)+wsgn*rcc*cos(abs(w));y1(end)+rcc*sin(abs(w))];
%center3 = [x1(end)-rcc*sin(pi/2-w);y1(end)+rcc*cos(pi/2-w)];
center3 = [x1(end)+h*cot(gamma+x)-rcc*cos(pi/2-gamma);y1(end)-rcc*sin(abs(w))];
%center3=[x1(end)-rcc*sin(3*pi/4-gamma);y1(end)+rcc*sin(3*pi/4-gamma)];
x3 = center3(1) + rcc*cos(th3);
y3 = center3(2) - rcc*sin(th3);
plot(x3,y3,'-g','linewidth',2);

plot(linspace(-4,6,100),y2(1)*ones(100,1),'-k','linewidth',2);
scatter(x1(1),y1(1),100,'bo','fill');
scatter(x2(1),y2(1),100,'bo','fill');
scatter(x1(end),y1(end),100,'bo','fill');
set(gca,'fontsize',20,'fontname','Times New Roman'); box on;
set(gcf,'color','w');
title({['T_{1} = ', num2str(t1), ' Alpha_{1} = ', num2str(a1*180/pi),' Beta_{1} = ', num2str(real(b1)*180/pi)],['T_{2} = ', num2str(t2), ' Alpha_{2} = ', num2str(a2*180/pi),' Beta_{2} = ',num2str(real(b2)*180/pi)],['Gamma = ', num2str(real(gamma)*180/pi),', w = ', num2str(real(w)*180/pi),', v = ', num2str(v),', h = ', num2str(h),', t_{cc} = ',num2str(tcc)]});
xlim([-4 4]); ylim([-0.5 1.5]);
pbaspect([8 2 1]);

% Compute area of trailer cell
xT = [x2,fliplr(x3),x2(1)];
yT = [y2,fliplr(y3),y2(1)];
aT = polyarea(xT,yT);

% Compute area of leader cell
xL = [x1,fliplr(x3),x1(1)];
yL = [y1,fliplr(y3),y1(1)];
aL = polyarea(xL,yL);

sprintf('Area of leader cell %.4f and trailer cell %.4f',aL,aT)


