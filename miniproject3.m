clear all; 
close all;

ax = 0.02;
ay = 0.01;
bx = 0.04;
bly = 0.01;
bhy = 0.015;
VA = 1.0;
VB = -1.0;
d = 0.001;
EPS0 = 8.8542*10^(-12);
presidue = 1e-4;

XAxis=bx;
YAxis=bly + ay +bhy;

NumX = round(XAxis/d) + 1;
NumY = round(YAxis/d) + 1;

distanceX = XAxis/(NumX-1);
distanceY = YAxis/(NumY-1);
x = linspace(0, XAxis, NumX);
y = linspace(0, YAxis, NumY);

left_gap = (bx-ax)/2;
xL = left_gap; 
xR = XAxis - left_gap; 
yB = bly; 
yT = bly + ay;

iLeft = round(xL/distanceX)+1;
iRight = round(xR/distanceX)+1;
jBottom = round(yB/distanceY)+1;
jTop = round(yT/distanceY)+1;

potential = zeros(NumX, NumY);

potential(1,:)=VB;
potential(end,:)=VB;
potential(:,1)=VB;
potential(:,end)=VB;

inner=false(NumX,NumY);
inner(iLeft:iRight,jBottom:jTop)=true;
potential(inner) = VA;

dielectric = true(NumX, NumY);
dielectric(1,:) = false;
dielectric(end,:) = false;
dielectric(:,1) = false;
dielectric(:,end) = false;
dielectric(inner) = false;

Maxiterations = 1e6;
iterations = 0;


while true
    iterations = iterations + 1;
    potentialNew = potential;
    for i = 2:NumX-1
        for j = 2:NumY-1
            if dielectric(i,j)
                potentialNew(i,j)=0.25*(potential(i+1,j)+potential(i-1,j)+potential(i,j+1)+potential(i,j-1));
            end
        end
    end
    delta = abs(potentialNew(dielectric)-potential(dielectric));
    Maxchange = max(delta);
    if (Maxchange<presidue)||(iterations>=Maxiterations)
        potential=potentialNew;
        break;
    end
    potential=potentialNew;
end

%Part A
[Xm, Ym] = meshgrid(x, y);
figure('Name','Part (A) 3D Surface of V');
S = surf(Xm, Ym, potential', 'EdgeColor','none', 'FaceColor','interp');
colormap(turbo(256));       
clim([VB VA]);             
shading interp;
xlabel('x [m]'); 
ylabel('y [m]'); 
zlabel('V [V]');
title('Part (A) 3D surface of V');
view(35,30); 
grid on;

%Part B
Ex=zeros(NumX, NumY);
Ey=zeros(NumX, NumY);
Ex(2:NumX-1,:)=-(potential(3:NumX,:)-potential(1:NumX-2,:))/(2*distanceX);
Ey(:,2:NumY-1)=-(potential(:,3:NumY)-potential(:,1:NumY-2))/(2*distanceY);
Ex(1,:)=-(potential(2,:)-potential(1,:))/distanceX;
Ex(end,:)=-(potential(end,:)-potential(end-1,:))/distanceX;
Ey(:,1)=-(potential(:,2)-potential(:,1))/distanceY;
Ey(:,end)=-(potential(:,end)-potential(:,end-1))/distanceY;
Ex(inner)=NaN;  
Ey(inner)=NaN;
figure('Name','Part (B) Electric Field');
quiver(Xm, Ym, Ex', Ey', 'AutoScale','on');
axis equal tight; 
grid on;
xlabel('x [m]'); ylabel('y [m]');
title('Part (B) Electric Field');

%Part C
jInnerTop=jTop;         
jAboveInner=jTop + 1;    
iSpanInner=iLeft:iRight; 
EinnerTop=-(potential(iSpanInner,jAboveInner)-potential(iSpanInner,jInnerTop))/distanceY;
sigmaInnerTop=EPS0*EinnerTop;
xInner=x(iSpanInner);
jOuterBottom=1;
jAboveOuter=2;
iSpanOuter=1:NumX;
EouterBottom=-(potential(iSpanOuter,jAboveOuter)-potential(iSpanOuter,jOuterBottom))/distanceY;
sigmaOuterBottom=EPS0*EouterBottom;
xOuter=x(iSpanOuter);
figure('Name','Part (C) Surface Charge Densities');
plot(xInner, sigmaInnerTop, 'LineWidth', 1.8); hold on;
plot(xOuter, sigmaOuterBottom, 'LineWidth', 1.8);
grid on; xlabel('x [m]'); ylabel('\sigma [C/m^2]');
title('Part (C) Surface Charge Density: Inner Top (Blue) vs Outer Bottom (Orange)');
xlim([0, XAxis]);
