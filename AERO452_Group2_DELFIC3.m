%% AERO 452 - PROJECT 2

%   Kate Parkinson 
%   Martin Kamme


%% HOUSE KEEPING

clear, clc, close all;

% Nomenclature
%   R = vector
%   r = magnitude
%   [] = units
%   R_X = vector in X frame



% Constants 
   global mue Re
     mue = 398600; % [kg^3/m^2]
     Re = 6378;    % [km]

 %% PART 1    
muearth = 398600; 
% Satellite Information
% TLE: 11/17
% DELFI-C3 (DO-64)        
TLE = [32789,  08021,     19320.66665100,  .00001377,  0,         .95251e-4,  0,                9999;
        32789,  97.4474,   5.6562,          .0011216,   335.1825,  24.8862,    15.072584486,     29605];

% Initial Orbit/Position
[COES0, JD0] = tle2coes(TLE);
ecc0 = COES0(1);
inc0 = deg2rad(COES0(2)); %
raan0 = deg2rad(COES0(3)); %
omega0 = deg2rad(COES0(4)); %
theta0 = deg2rad(COES0(5)); %     
a0 = COES0(6);
h0 = COES0(7);
T0 = ((2*pi)/sqrt(muearth))*a0^(3/2); % Period of the original orbit

% Geometry/Mass
mDelfi = 3; % kg
rDelfi = (34/2)*(1/100); % m

%Propogating forward with VoP
timelong = 365*24*60*60; % seconds
timeshort = 30*24*60*60; % seconds

tspan = linspace(0,timelong,400);
tspanshort = linspace(0,timeshort,400);

coesin = [h0;
          ecc0;
          theta0;
          raan0;
          inc0;
          omega0];
      

options = odeset('RelTol',1e-8,'AbsTol', 1e-8,'events',@Terminate);


[tnewshort,coesoutshort] = ode45(@VoP,tspanshort,coesin,options,JD0);

[tnew,coesout] = ode45(@VoP,tspan,coesin,options,JD0);

    



%% Plot VoP 

% Short term---------------------------------------------------------------
timeplotshort = tnewshort/(24*60*60);
hplotshort = coesoutshort(:,1);
eccplotshort = smooth((coesoutshort(:,2)));
thetaplotshort = rad2deg(coesoutshort(:,3));
RAANplotshort = rad2deg(coesoutshort(:,4));
incplotshort = rad2deg(coesoutshort(:,5));
perplotshort = rad2deg(coesoutshort(:,6));

for i = 1:length(tnewshort)
    placeholder = perplotshort(i);
    while placeholder >= 360
        placeholder = placeholder - 360;
    end
    perplotshort(i) = placeholder;
end




for i = 1:length(timeplotshort)
    rashort(i) = ((((hplotshort(i))^2)/mue)/(1+eccplotshort(i)*cosd(180))) - Re;
    rpshort(i) = ((((hplotshort(i))^2)/mue)/(1+eccplotshort(i)*cosd(0))) - Re;
end


figure
title('Delfi-C3 Short Term')
hold on
plot(timeplotshort,rashort,'b','linewidth',2)
plot(timeplotshort,rpshort,'r','linewidth',2)
grid on 
grid minor
xlabel('Time (days)')
ylabel('Altitude (km)')
legend('ra','rp')
hold off

figure
subplot(4,1,1)
plot(timeplotshort,eccplotshort)
grid on
grid minor
xlabel('Time (days)')
ylabel('ecc')
title('Delfi-C3 Short Term')

 
subplot(4,1,2)
plot(timeplotshort,incplotshort)
grid on
grid minor
xlabel('Time (days)')
ylabel(['Inc ' char(176)])

subplot(4,1,3)
plot(timeplotshort,RAANplotshort)
grid on
grid minor
xlabel('Time (days)')
ylabel(['RAAN ' char(176)])

subplot(4,1,4)
plot(timeplotshort,perplotshort)
grid on
grid minor
xlabel('Time (days)')
ylabel(['\omega' char(176)])


% Long Term----------------------------------------------------------------
timeplot = tnew/(24*60*60);
hplot = coesout(:,1);
eccplot = smooth((coesout(:,2)));
thetaplot = rad2deg(coesout(:,3));
RAANplot = rad2deg(coesout(:,4));
incplot = rad2deg(coesout(:,5));
perplot = rad2deg(coesout(:,6));


for i = 1:length(tnewshort)
    placeholder = perplot(i);
    while placeholder >= 360
        placeholder = placeholder - 360;
    end
    perplot(i) = placeholder;
end

for i = 1:length(timeplot)
    ra(i) = ((((hplot(i))^2)/mue)/(1+eccplot(i)*cosd(180))) - Re;
    rp(i) = ((((hplot(i))^2)/mue)/(1+eccplot(i)*cosd(0))) - Re;
end


figure
title('Delfi-C3 Long Term')
hold on
plot(timeplot,ra,'b','linewidth',2)
plot(timeplot,rp,'r','linewidth',2)
grid on 
grid minor
xlabel('Time (days)')
ylabel('Altitude (km)')
legend('ra','rp')
hold off

figure
subplot(4,1,1)
plot(timeplot,eccplot)
grid on
grid minor
xlabel('Time (days)')
ylabel('ecc')
title('Delfi-C3 Long Term')

 
subplot(4,1,2)
plot(timeplot,incplot)
grid on
grid minor
xlabel('Time (days)')
ylabel(['Inc ' char(176)])

subplot(4,1,3)
plot(timeplot,RAANplot)
grid on
grid minor
xlabel('Time (days)')
ylabel(['RAAN ' char(176)])

subplot(4,1,4)
plot(timeplot,perplot)
grid on
grid minor
xlabel('Time (days)')
ylabel(['\omega' char(176)])



%% Comparison with osculating orbit
close all
clear all
clc

global mue Re
     mue = 398600; % [kg^3/m^2]
     Re = 6378;    % [km]
     
muearth = 398600; % km^3/s^2

% Satellite Information
% TLE: 11/17
% DELFI-C3 (DO-64)        
TLE = [32789,  08021,     19320.66665100,  .00001377,  0,         .95251e-4,  0,                9999;
       32789,  97.4474,   5.6562,          .0011216,   335.1825,  24.8862,    15.072584486,     29605];

% Initial Orbit/Position
[COES0, JD0] = tle2coes(TLE);
ecc0 = COES0(1);
inc0 = deg2rad(COES0(2));
raan0 = deg2rad(COES0(3));
omega0 = deg2rad(COES0(4));
theta0 = deg2rad(COES0(5));
a0 = COES0(6);
h0 = COES0(7);
T0 = ((2*pi)/sqrt(muearth))*a0^(3/2); % Period of the original orbit

[R,V] = coes2rvKAMME(h0,inc0,raan0,ecc0,omega0,theta0,muearth);


% Initial Conditions for While Loop
timeFINISH = 7*12*24*60*60; % 12 weeks 
dvTsum = 0;
coesin = [h0;
          ecc0;
          theta0;
          raan0;
          inc0;
          omega0];
state = [R;V];
t = 0;
options = odeset('RelTol',1e-8,'AbsTol', 1e-8);

[x,y,z] = sphere;
x = x*6378;
y = y*6378;
z = z*6378;
figure
hold on
surf(x,y,z)
axis equal
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

while t < timeFINISH
    
% Geometry/Mass
mDelfi = 3; % kg
rDelfi = (34/2)*(1/100); % m

% Correct every week
tprop = 7*24*60*60; % seconds 
tspan = linspace(0,tprop);

% Burn time for lambert
timeLAM = T0/2; 
tspanLAM = linspace(0,timeLAM);      


% VoP-----------
[tnew,coesout] = ode45(@VoP,tspan,coesin,options,JD0);

hf = coesout(length(coesout),1);
eccf = coesout(length(coesout),2);
thetaf = coesout(length(coesout),3);
RAANf = coesout(length(coesout),4);
incf = coesout(length(coesout),5);
perf = coesout(length(coesout),6);
[Rf,Vf] = coes2rvKAMME(hf,incf,RAANf,eccf,perf,thetaf,muearth);


% 2body----------
[tosc,stateoutosc] = ode45(@twoBody,tspan,state,options);

% Before Lamberts burn
Rosc = [stateoutosc(length(stateoutosc),1:3)]';
Vosc = [stateoutosc(length(stateoutosc),4:6)]';
state = [Rosc;Vosc];


% After Lamberts Burn
[toscf,stateoutoscf] = ode45(@twoBody,tspanLAM,state,options);
Roscf = [stateoutoscf(length(stateoutoscf),1:3)]';
Voscf = [stateoutoscf(length(stateoutoscf),4:6)]';

% Lamberts Burn
string = 'retro';
[V1, V2] = lambert(Rf, Roscf, timeLAM, string);


% DeltaV 
dvleave = abs(norm(V1 - Vf));
dvarrive = abs(norm(Voscf - V2));
dvT = dvleave + dvarrive;

% Plot
stateplot = [Rosc;V1];
[toutplot,stateoutplot] = ode45(@twoBody,tspanLAM,stateplot,options);
plot3(stateoutplot(1,1),stateoutplot(1,2),stateoutplot(1,3),'r*','MarkerSize',10)
plot3(stateoutplot(length(stateoutplot),1),stateoutplot(length(stateoutplot),2),stateoutplot(length(stateoutplot),3),'g*','MarkerSize',10)
plot3(stateoutplot(:,1),stateoutplot(:,2),stateoutplot(:,3))
% plot3(stateoutoscf(:,1),stateoutoscf(:,2),stateoutoscf(:,3))
% plot3(stateoutplotactual(:,1),stateoutplotactual(:,2),stateoutplotactual(:,3))

% Conditions for beginning of next iteration
state = [Roscf;Voscf];
[a0,ecc0,inc0,raan0,omega0,ME,truLon,argLat,lonPer,p] = rv2orb(Roscf,Voscf,muearth);
h0 = sqrt(muearth*a0*(1-ecc0^2));

% Calculate the eccentric anomaly using Mean anomaly
err = 1e-8;            %Calculation Error
E0 = ME;
tcalc =1;
itt = 0;
while(tcalc) 
    E =  ME + ecc0*sin(E0);
    if ( abs(E - E0) < err)
        tcalc = 0;
    end
    E0 = E;
    itt = itt+1;
end

% Calc true anomaly with newton's iterative scheme 
rat = 25;
count = 0;
Eanom = E;
Manom = ME;
while abs(rat) > 10^-8 && count < 10000
    Eanom = Eanom - (-Manom + Eanom - ecc0*sin(Eanom))/(1-ecc0*cos(Eanom));
    count = count +1;
    rat = (-Manom + Eanom - ecc0*sin(Eanom))/(1-ecc0*cos(Eanom));
end

theta0 = 2*atan((tan(Eanom/2))*sqrt((1+ecc0)/(1-ecc0))); 

coesin = [h0;
          ecc0;
          theta0;
          raan0;
          inc0;
          omega0];
    
t = t + tprop + timeLAM;
dvTsum = dvTsum + dvT;

end
legend('Earth','Leave','Arrive')


%% ~~~~~~~ FUNCTIONS ~~~~~~~ %%

%% TLE to RV

function[R,V, Jdate, P, n, COES] = TLE_to_RV(TLE)
    % Inputs: TLE of target
    % Outputs: R,V of target in ECI

global mue
    
% Take TLE and separate into classic orbital elements and Julian date
% Format: COES =[Semi-Mjr Axis, ecc, inc(deg), RAAN(deg) w(deg), true anomaly(deg)];
[COEs, Jdate] = tle2coes(TLE);
    

ecc = COEs(1);
inc = rad2deg(COEs(2));
RAAN = rad2deg(COEs(3));
arg_p = rad2deg(COEs(4));
tru_a = rad2deg(COEs(5));
a = COEs(6);
h = sqrt(a*mue*(1-ecc^2));
 
% Period of target
P = (2*pi*a^1.5)/sqrt(mue);
 
% Mean Motion of target
n = sqrt(mue/a^3);
    
[R,V] = coes2rv(ecc, inc, RAAN, arg_p, tru_a, h);    

COES = [ecc, inc, RAAN, arg_p, tru_a, a, h];

end

function [COEs, date1] = tle2coes(TLE)

global mue

epoch =  TLE(1,3)*24*3600;   % Epoch Date and Julian Date Fraction
Db    = TLE(1,6);            % Ballistic Coefficient
inc   = TLE(2,2);            % Inclination [deg]
RAAN  = TLE(2,3);            % Right Ascension of the Ascending Node [deg] 
ecc     = TLE(2,4);          % Eccentricity 
w     = TLE(2,5);            % Argument of periapsis [deg]
M     = TLE(2,6);            % Mean anomaly [deg]
n     = TLE(2,7);            % Mean motion [Revs per day]


% Orbital elements
SMA = (mue/(n*2*pi/(24*3600))^2)^(1/3);     % Semi-major axis [km]    
h = sqrt(SMA*mue*(1-ecc^2));

% Calculate the eccentric anomaly using Mean anomaly
err = 1e-8;            %Calculation Error
E0 = M;
t =1;
itt = 0;
while(t) 
    E =  M + ecc*sind(E0);
    if ( abs(E - E0) < err)
        t = 0;
    end
    E0 = E;
    itt = itt+1;
end

% Calc true anomaly with newton's iterative scheme 
rat = 25;
count = 0;
Eanom = E;
Manom = M;
while abs(rat) > 10^-8 && count < 10000
    Eanom = Eanom - (-Manom + Eanom - ecc*sind(Eanom))/(1-ecc*cosd(Eanom));
    count = count +1;
    rat = (-Manom + Eanom - ecc*sind(Eanom))/(1-ecc*cosd(Eanom));
end
true_a = 2*atand((tand(Eanom/2))*sqrt((1+ecc)/(1-ecc))); 
    
COEs = [ecc, inc, RAAN, w, true_a, SMA,h];
date1 = TLE(1,3);
days = TLE(1,3) - 19000;  % Epoch Date and Julian Date Fraction
yr = 2450000;
date1 = yr + 8485 + days;
end


%% Standard Functions

function [r_ECI,v_ECI] = coes2rv(ecc, incl, RAAN, arg_p, tru_a, h)
    %INPUTs: COES in deg, h is either magnitude OR vector
    %OUTPUTs: R, V in ECI

global mue   
h = norm(h);
    
% Finding R and V vectors in perifocal ref frame from COES
r_mag = ((norm(h)^2)/mue)*(1/(1+ecc*cos(tru_a)));

R_PERI = r_mag*[cos(tru_a);sin(tru_a); 0];
V_PERI = (mue/norm(h))*[-sin(tru_a); (ecc + cos(tru_a));0];    

% Rotation matrices for each angle
R3_w = [cos(arg_p) sin(arg_p) 0; -sin(arg_p) cos(arg_p) 0; 0 0 1];
R1_inc = [1 0 0; 0 cos(incl) sin(incl); 0 -sin(incl) cos(incl)];
R3_raan = [cos(RAAN) sin(RAAN) 0; -sin(RAAN) cos(RAAN) 0; 0 0 1];

% Rotation matrix Q(peri<-eci) using 3-1-3 rotaion sequence
Q = R3_w*R1_inc*R3_raan;

% Find r and v vectors in ECI using Q' to go from Peri->ECI
r_ECI = Q'*R_PERI;
v_ECI = Q'*V_PERI;

end

function [R,V] = coes2rvKAMME(h,inc,RAAN,e,per,theta,muearth)
% h [km^2/s] Specific angular momentum
% i [rad] Inclination
% RAAN [rad] Right ascension (RA) of the ascending node
% e Eccentricity
% per [rad] Argument of perigee
% theta [rad] True anomaly
% muearth = 398600; Earth’s gravitational parameter [km^3/s^2]

% State Vectors in Perifocal coordinates
rx = h^2/muearth*(1/(1 + e*cos(theta)))*[cos(theta);sin(theta);0];
vx = muearth/h*[-sin(theta); (e +cos(theta));0];

% Direction cosine matrix
DCM = [cos(per), sin(per),0;-sin(per),cos(per),0;0,0,1]*...
 [1,0,0;0,cos(inc),sin(inc);0,-sin(inc),cos(inc)]*...
 [cos(RAAN), sin(RAAN),0;-sin(RAAN),cos(RAAN),0;0,0,1];

% Transformation Matrix
Dcm = inv(DCM);

% ECI R
R = Dcm*rx;

% ECI V
V = Dcm*vx;

end

function [a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(r,v,mu)
%----------------------- Begin Code Sequence -----------------------------%
% Purpose:                                                                %
% Convert a given set of state vectors in ECI reference frame to orbital  %
% elements.                                                               %
%-------------------------------------------------------------------------%
%                                                                         %
% Inputs:                                                                 %
%--------                                                                  
%r_ECI                  [3 x N]                         Position Vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
%v_ECI                  [3 x N]                         Velocity vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
%mu                     double                          Gravitational Constant
%                                                       Defaults to Earth if
%                                                       not specified
% Outputs:
%---------                                                                %
%a                      [1 x N]                         Semi-Major Axis
%                                                       (km)
%
%eMag                   [1 x N]                         Eccentricity
%                                                       (unitless)
%
%i                      [1 x N]                         inclination
%                                                       (radians)
%
%O                      [1 x N]                         Right Ascention of
%                                                       the ascending node
%                                                       (radians)
%
%o                      [1 x N]                         Argument of perigee
%                                                       (radians)
%
%M                      [1 x N]                         Mean Anomaly
%                                                       (radians)
%
%truLon                 [1 x N]                         True Longitude
%                                                       (radians)
%
%argLat                 [1 x N]                         Argument of Latitude
%                                                       (radians)
%
%lonPer                 [1 x N]                         Longitude of Periapse
%                                                       (radians)
%
%p                      [1 x N]                         Semilatus Rectum
%                                                       (km)
%
% References:
%-------------
%Vallado,D. Fundamentals of Astrodynamics and Applications. 2007.
%
% Function Dependencies:
%------------------
%None
%------------------------------------------------------------------       %
% Programed by Darin Koblick  03-04-2012                                  %
% Updated to address circular equatorial orbits       12/12/2013          %
%------------------------------------------------------------------       %

if ~exist('mu','var');  t = getConst(); mu = t.Earth.Mu; end
%Specific angular momentum
h = cross(r,v);
n = cross(repmat([0;0;1],[1,size(r,2)]),h); nMag = sqrt(sum(n.^2,1));
vMag = sqrt(sum(v.^2,1)); 
rMag = sqrt(sum(r.^2,1)); 
hMag = sqrt(sum(h.^2,1));
e = (1./mu).*(bsxfun(@times,(vMag.^2 - mu./rMag),r) - bsxfun(@times,dot(r,v),v)); 
eMag = sqrt(sum(e.^2,1));
zeta = (vMag.^2)./2 - mu./rMag;

%Special Procedure when we have a parabolic orbit
idx = eMag ~= 1;
a = NaN(size(eMag));
p = NaN(size(eMag));
if any(idx)
    a(idx) = -mu./(2.*zeta(idx)); 
    p = a(idx).*(1-eMag(idx).^2); 
else
    a(idx) = Inf; 
    p(idx) = (hMag(idx).^2)./mu; 
end

%Compute the angles
i = acos(h(3,:)./hMag); 
O = acos(n(1,:)./nMag);
o = acos(dot(n,e)./(nMag.*eMag));
nu = acos(dot(e,r)./(eMag.*rMag));
lonPer = acos(e(1,:)./eMag);
argLat = acos(dot(n,r)./(nMag.*rMag));
truLon = acos(r(1,:)./rMag);

%Account for those cases where satellite is in circular orbit
         O(n(1,:) == 0) = 0;
       o(dot(n,e) == 0) = 0;
    lonPer(e(1,:) == 0) = 0;
      nu(dot(e,r) == 0) = 0;
  argLat(dot(n,r) == 0) = 0;
  
%Apply Quadrant Checks to All Determined Angles
idx = n(2,:) < 0; if any(idx);  O(idx) = 2*pi - O(idx);  end
idx = e(3,:) < 0; if any(idx); o(idx) = 2*pi - o(idx); end
idx = dot(r,v) < 0; if any(idx); nu(idx) = 2*pi - nu(idx); end
idx = e(2,:) < 0; if any(idx); lonPer(idx) = 2*pi-lonPer(idx);  end
idx = r(3,:) < 0; if any(idx); argLat(idx) = 2*pi - argLat(idx); end
idx = r(2,:) < 0; if any(idx); truLon(idx) = 2*pi - truLon(idx); end
end

%% Propagation 

% VoP
function [coesdot] = VoP(t,coesin,JD0)
% Coesin includes the following
% coesin(1) = h
% coesin(2) = ecc
% coesin(3) = true anomaly
% coesin(4) = RAAN
% coesin(5) = inc
% coesin(6) = per
A = pi*((34/2)*(1/100))^2; % m
m = 3; % kg
Cr = 1.2;

hours = 3600; %Hours to seconds
days = 24*hours; %Days to seconds

muearth = 398600;
mu3 = 4903;
rearth = 6378;

J2 = .00108263;
J3 = (-2.33936*10^-3)*J2;
J4 = (-1.49601*10^-3)*J2;
J5 = (-.20995*10^-3)*J2;
J6 = (.49941*10^-3)*J2;


h = coesin(1);
ecc = coesin(2);
theta = coesin(3);
RAAN = coesin(4);
inc = coesin(5);
per = coesin(6);

[R,V] = coes2rvKAMME(h,inc,RAAN,ecc,per,theta,muearth);
x = R(1);
y = R(2);
z = R(3);

r = norm(R);

%...Obtain the unit vectors of the rsw system:
r = norm(R);
ur = R/r; %radial
H = cross(R,V);
uh = H/norm(H); %normal
s = cross(uh, ur);
us = s/norm(s); %transverse

Qxr = [-sin(RAAN)*cos(inc)*sin(per+theta)+cos(RAAN)*cos(per+theta) , cos(RAAN)*cos(inc)*sin(per+theta)+sin(RAAN)*cos(per+theta), sin(inc)*sin(per+theta);
       -sin(RAAN)*cos(inc)*cos(per+theta)-cos(RAAN)*sin(per+theta) , cos(RAAN)*cos(inc)*cos(per+theta)-sin(RAAN)*sin(per+theta), sin(inc)*cos(per+theta);
       sin(RAAN)*sin(inc) , -cos(RAAN)*sin(inc) , cos(inc)];
  
% J3   
aI3 = -((5*J3*muearth*(rearth^3)*x)/(2*r^7))*((3*z) - ((7*(z^3))/(r^2)));
aJ3 = -((5*J3*muearth*(rearth^3)*y)/(2*r^7))*((3*z) - ((7*(z^3))/(r^2)));
aK3 = -((5*J3*muearth*(rearth^3))/(2*r^7))*((6*(z^2)) - ((7*(z^4))/(r^2))-(3/5)*(r^2));

P3 = [aI3;
     aJ3;
     aK3];
 
Prsw3  = Qxr*P3;

Pr3 = Prsw3(1);
Ps3 = Prsw3(2);
Pw3 = Prsw3(3);

% J4 
aI4 = ((15*J4*muearth*(rearth^4)*x)/(8*r^7))*(1-((14*z^2)/(r^2)) + ((21*z^4)/(r^4)));
aJ4 = ((15*J4*muearth*(rearth^4)*y)/(8*r^7))*(1-((14*z^2)/(r^2)) + ((21*z^4)/(r^4)));
aK4 = ((15*J4*muearth*(rearth^4)*z)/(8*r^7))*(5-((70*z^2)/(3*(r^2))) + ((21*z^4)/(r^4)));

P4 = [aI4;
     aJ4;
     aK4];
 
Prsw4  = Qxr*P4;

Pr4 = Prsw4(1);
Ps4 = Prsw4(2);
Pw4 = Prsw4(3);


% J5 
aI5 = ((3*J5*muearth*(rearth^5)*x*z)/(8*(r^9)))*(35 - (210*((z^2)/(r^2))) + (231*((z^4)/(r^4))));
aJ5 = ((3*J5*muearth*(rearth^5)*y*z)/(8*(r^9)))*(35 - (210*((z^2)/(r^2))) + (231*((z^4)/(r^4))));
aK5 = ((3*J5*muearth*(rearth^5)*z*z)/(8*(r^9)))*(105 - (315*((z^2)/(r^2))) + (231*((z^4)/(r^4)))) - ((15*J5*muearth*(rearth^5))/(8*r^7));


P5 = [aI5;
     aJ5;
     aK5];
 
Prsw5  = Qxr*P5;

Pr5 = Prsw5(1);
Ps5 = Prsw5(2);
Pw5 = Prsw5(3);

% J6
aI6 = -((J6*muearth*(rearth^6)*x)/(16*(r^9)))*(35 - (945*((z^2)/(r^2))) + (3465*((z^4)/(r^4))) - (3003*((z^6)/(r^6))));
aJ6 = -((J6*muearth*(rearth^6)*y)/(16*(r^9)))*(35 - (945*((z^2)/(r^2))) + (3465*((z^4)/(r^4))) - (3003*((z^6)/(r^6))));
aK6 = -((J6*muearth*(rearth^6)*z)/(16*(r^9)))*(245 - (2205*((z^2)/(r^2))) + (4851*((z^4)/(r^4))) - (3003*((z^6)/(r^6))));

P6 = [aI6;
     aJ6;
     aK6];
 
Prsw6  = Qxr*P6;

Pr6 = Prsw6(1);
Ps6 = Prsw6(2);
Pw6 = Prsw6(3);

% Moon Perturbation

%...Update the Julian Day:
JD = JD0 + t/days;

%...Find and normalize the position vector of the moon:
R_m = lunar_position(JD);
R_m = R_m'; % convert to columns
r_m = norm(R_m);
R_rel = R_m - R; %R_rel = position vector of moon wrt satellite
r_rel = norm(R_rel);

% Some garbage Curtis Includes for whatever reason
q = dot(R,(2*R_m - R))/r_m^2;
F = (q^2 - 3*q + 3)*q/(1 + (1-q)^1.5);

% Gravitational perturbation of the moon (Equation 12.117):
apMOON = mu3/r_rel^3*(F*R_m - R);

apMOONR = dot(apMOON,ur);
apMOONS = dot(apMOON,us);
apMOONW = dot(apMOON,uh);

% SRP 
[Asrp] = SRP(R,A,m,Cr,JD);
AsrpRSW = Qxr*Asrp;

% drag
ap = dragKamme(R,V);
apRSW = Qxr*(ap/1000);

% All perturbations in RSW
Pr = apMOONR;% -(3/2)*((J2*muearth*(rearth^2))/(r^4))*(1-3*((sin(inc))^2)*((sin(per+theta))^2)) + Pr3 + Pr4 + Pr5 + Pr6 + AsrpRSW(1) + apRSW(1); %

Ps = apMOONS;% -(3/2)*((J2*muearth*(rearth^2))/(r^4))*((sin(inc))^2)*sin(2*(per+theta)) + Ps3 + Ps4 + Ps5 + Ps6  + AsrpRSW(2) + apRSW(2); % 

Pw = apMOONW; % -(3/2)*((J2*muearth*(rearth^2))/(r^4))*(sin(2*(inc))*sin(per+theta)) + Pw3 + Pw4 + Pw5 + Pw6  + AsrpRSW(3) + apRSW(3); %


% Coesdot 
    dh = r*Ps;
    
    decc = (h/muearth)*sin(theta)*Pr + (1/(muearth*h))*((h^2 + muearth*r)*cos(theta) + muearth*ecc*r)*Ps;
    
    dtheta = (h/r^2) + (1/(ecc*h))*(((h^2)/muearth)*Pr*cos(theta) - (((h^2)/muearth)+r)*Ps*sin(theta));
    
    dRAAN = ((r*sin(per+theta))/(h*sin(inc)))*Pw;
    
    dinc = (r/h)*cos(per+theta)*Pw;
    
    dper = -((r*sin(per+theta))/(h*tan(inc)))*Pw - (1/(ecc*h))*(((h^2)/muearth)*Pr*cos(theta) - (((h^2)/muearth)+r)*Ps*sin(theta));
    
    
coesdot = [dh;
           decc;
           dtheta;
           dRAAN;
           dinc;
           dper];
       
end

% 2 Body
function [statedot] = twoBody(t,state)
% state includes the following: 
% state(1) = rAx eci
% state(2) = rAy eci
% state(3) = rAz eci
% state(4) = vAx eci
% state(5) = vAy eci
% state(6) = vAz eci
muearth = 398600;

R = state(1:3);
V = state(4:6);

r = norm(R);

a0 = -muearth*R/r^3;
a = a0;
 

statedot =  [V;a];  % Propagation of R and V vectors

end


%% Perturbations
function [Asrp] = SRP(R,A,m,Cr,JD)
  
    % A = Area facing the sun
    % Cr = reflectivity of s/c
    % Psr = Solar Radiation Pressure
    % m =  S/c mass
    % Rsat_sun = vector from satellite to sun
    % F = shadown function (0,1)
    
    R_sat = R;
    
    Psr = (1370/(3e8))/1000; %E/c [Pa]
       
    [~,~, R_sun] = solar_position(JD);
    
    F = los(R_sat, R_sun);
    
    rsun = R_sun/norm(R_sun);

    Asrp = -(Psr*Cr*A/m)*rsun*F;
end

function [adrag] = dragKamme(R,V)
Cd = 2.2;

% Geometry/Mass
mDelfi = 3; % kg
rDelfi = (34/2)*(1/100); % m

d = 1; %m
A = pi*(rDelfi)^2; %m^2
m = mDelfi; %kg
omegaEarth = [0;0;72.9211*10^-6]; % rad/s

height = (norm(R) - 6378); %km 
rho = atmosphere(height);
vrel = V - cross(omegaEarth,R);

adrag = (-1/2)*((Cd*A)/m)*rho*((1000*norm(vrel))^2)*(vrel/norm(vrel));
end

% Exponential Density Model (Cite: Curtis)
function [rho] = atmosphere(z)

%...Geometric altitudes (km):
h = ...
[ 0 25 30 40 50 60 70 ...
80 90 100 110 120 130 140 ...
150 180 200 250 300 350 400 ...
450 500 600 700 800 900 1000];

%...Corresponding densities (kg/m^3) from USSA76:
r = ...
[1.225 4.008e-2 1.841e-2 3.996e-3 1.027e-3 3.097e-4 8.283e-5 ...
1.846e-5 3.416e-6 5.606e-7 9.708e-8 2.222e-8 8.152e-9 3.831e-9 ...
2.076e-9 5.194e-10 2.541e-10 6.073e-11 1.916e-11 7.014e-12 2.803e-12 ...
1.184e-12 5.215e-13 1.137e-13 3.070e-14 1.136e-14 5.759e-15 3.561e-15];

%...Scale heights (km):
H = ...
[ 7.310 6.427 6.546 7.360 8.342 7.583 6.661 ...
5.927 5.533 5.703 6.782 9.973 13.243 16.322 ...
21.652 27.974 34.934 43.342 49.755 54.513 58.019 ...
60.980 65.654 76.377 100.587 147.203 208.020];

%...Handle altitudes outside of the range:
if z > 1000
    z = 1000;
elseif z < 0
    z = 0;
end

%...Determine the interpolation interval:
for j = 1:27
    if z >= h(j) && z < h(j+1)
        i = j;
    end

end

if z == 1000
    i = 27;
end
%...Exponential interpolation:
rho = r(i)*exp(-(z - h(i))/H(i));
end


%% Useful external Functions

% Moon Position (Cite: Curtis)
function r_moon = lunar_position(jd)
%Calculates the geocentric equatorial position vector of the moon given the Julian day.

RE = 6378;

%...Time in centuries since J2000:
T = (jd - 2451545)/36525;

%...Ecliptic longitude (deg):
e_long = 218.32 + 481267.881*T ...
+ 6.29*sind(135.0 + 477198.87*T) - 1.27*sind(259.3 - 413335.36*T)...
+ 0.66*sind(235.7 + 890534.22*T) + 0.21*sind(269.9 + 954397.74*T)...
- 0.19*sind(357.5 + 35999.05*T) - 0.11*sind(186.5 + 966404.03*T);
e_long = mod(e_long,360);


%...Ecliptic latitude (deg):
e_lat = 5.13*sind( 93.3 + 483202.02*T) + 0.28*sind(228.2 + 960400.89*T)...
- 0.28*sind(318.3 + 6003.15*T) - 0.17*sind(217.6 - 407332.21*T);
e_lat = mod(e_lat,360);


%...Horizontal parallax (deg):
h_par = 0.9508 ...
+ 0.0518*cosd(135.0 + 477198.87*T) + 0.0095*cosd(259.3 - 413335.36*T)...
+ 0.0078*cosd(235.7 + 890534.22*T) + 0.0028*cosd(269.9 + 954397.74*T);
h_par = mod(h_par,360);


%...Angle between earth’s orbit and its equator (deg):
obliquity = 23.439291 - 0.0130042*T;


%...Direction cosines of the moon’s geocentric equatorial position vector:
l = cosd(e_lat) * cosd(e_long);
m = cosd(obliquity)*cosd(e_lat)*sind(e_long) - sind(obliquity)*sind(e_lat);
n = sind(obliquity)*cosd(e_lat)*sind(e_long) + cosd(obliquity)*sind(e_lat);

%...Direction cosines of the moon’s geocentric equatorial position vector:
l = cosd(e_lat) * cosd(e_long);
m = cosd(obliquity)*cosd(e_lat)*sind(e_long) - sind(obliquity)*sind(e_lat);
n = sind(obliquity)*cosd(e_lat)*sind(e_long) + cosd(obliquity)*sind(e_lat);


%...Earth-moon distance (km):
dist = RE/sind(h_par);


%...Moon’s geocentric equatorial position vector (km):
r_moon = dist*[l m n];


end %lunar_position

% Eclipse (Cite: Curtis)
function light_switch = los(r_sat, r_sun)
%%CREDIT: CURTIS

% This function uses the ECI position vectors of the satellite (r_sat)
% and the sun (r_sun) to determine whether the earth is in the line of
% sight between the two.
%
% User M-functions required: None
%--------------------------------------------------------------------------
RE        = 6378;          %Earth's radius (km)
rsat      = norm(r_sat);
rsun      = norm(r_sun);

%...Angle between sun and satellite position vectore:  
theta     = acosd(dot(r_sat, r_sun)/rsat/rsun);

%...Angle between the satellite position vector and the radial to the point
%   of tangency with the earth of a line from the satellite:                                              
theta_sat = acosd(RE/rsat);

%...Angle between the sun position vector and the radial to the point
%   of tangency with the earth of a line from the sun: 
theta_sun = acosd(RE/rsun);

%...Determine whether a line from the sun to the satellite 
%   intersects the earth:
if theta_sat + theta_sun <= theta
    light_switch = 0;   %yes
else
    light_switch = 1;   %no
end

end %los

% Sun Location (Cite: Curtis)
function [lamda eps r_S] = solar_position(jd)
%CREDIT: CURTIS

    % This function alculates the geocentric equatorial position vector
    % of the sun, given the julian date.
    %
    % -------------------------------------------------------------------------
    %...Astronomical unit (km):
    AU    = 149597870.691;

    %...Julian days since J2000:
    n     = jd - 2451545;

    %...Julian centuries since J2000:
    cy    = n/36525;

    %...Mean anomaly (deg{:
    M     = 357.528 + 0.9856003*n;
    M     = mod(M,360);

    %...Mean longitude (deg):
    L     = 280.460 + 0.98564736*n;
    L     = mod(L,360);

    %...Apparent ecliptic longitude (deg):
    lamda = L + 1.915*sind(M) + 0.020*sind(2*M);
    lamda = mod(lamda,360);

    %...Obliquity of the ecliptic (deg):
    eps   = 23.439 - 0.0000004*n;

    %...Unit vector from earth to sun:
    u     = [cosd(lamda); sind(lamda)*cosd(eps); sind(lamda)*sind(eps)];

    %...Distance from earth to sun (km):
    rS    = (1.00014 - 0.01671*cosd(M) - 0.000140*cosd(2*M))*AU;

    %...Geocentric position vector (km):
    r_S   = rS*u;
    end %solar_position
 
% Lamberts
function [V1, V2] = lambert(R1, R2, t, string)

global mue

%...Magnitudes of R1 and R2:
r1 = norm(R1);
r2 = norm(R2);
c12 = cross(R1, R2);
theta = acos(dot(R1,R2)/r1/r2);

%...Determine whether the orbit is prograde or retrograde:
if nargin < 4 || (~strcmp(string,'retro') & (~strcmp(string,'pro')))
    string = 'pro';
    fprintf('\n ** Prograde trajectory assumed.\n')
end

if strcmp(string,'pro')
    if c12(3) <= 0
        theta = 2*pi - theta;
    end
elseif strcmp(string,'retro')
    if c12(3) >= 0
        theta = 2*pi - theta;
    end
end

%...Equation 5.35:
A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));

%...Determine approximately where F(z,t) changes sign, and
%...use that value of z as the starting value for Equation 5.45:
z = -100;
while F(z,t) < 0
    z = z + 0.1;
end

%...Set an error tolerance and a limit on the number of iterations:
tol = 1.e-8;
nmax = 5000;

%...Iterate on Equation 5.45 until z is determined to within the
%...error tolerance:
ratio = 1;
n = 0;
while (abs(ratio) > tol) & (n <= nmax)
    n = n + 1;
    ratio = F(z,t)/dFdz(z);
    z = z - ratio;
end

%...Report if the maximum number of iterations is exceeded:
if n >= nmax
fprintf('\n\n **Number of iterations exceeds %g \n\n ',nmax)
end

%...Equation 5.46a:
f = 1 - y(z)/r1;

%...Equation 5.46b:
g = A*sqrt(y(z)/mue);

%...Equation 5.46d:
gdot = 1 - y(z)/r2;

%...Equation 5.28:
V1 = 1/g*(R2 - f*R1);

%...Equation 5.29:
V2 = 1/g*(gdot*R2 - R1);

return

%...Equation 5.38:
function dum = y(z)
    dum = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
end

%...Equation 5.40:
function dum = F(z,t)
    dum = (y(z)/C(z))^1.5*S(z) + A*sqrt(y(z)) - sqrt(mue)*t;
end

%...Equation 5.43:
function dum = dFdz(z)
    if z == 0
        dum = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
    else
        dum = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) ...
               + 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z)) ...
               + A*sqrt(C(z)/y(z)));
    end
end

%...Stumpff functions:
function dum = C(z)
    dum = stumpC(z);
end

function dum = S(z)
    dum = stumpS(z);
end

end %lambert

function c = stumpC(z)

if z > 0
    c = (1 - cos(sqrt(z)))/z;
elseif z < 0
    c = (cosh(sqrt(-z)) - 1)/(-z);
else
c = 1/2;
end

end

function s = stumpS(z)
if z > 0
    s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
elseif z < 0
    s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
else
    s = 1/6;
end
end




%% Terminate

% Terminate ODE VOP
function [lookfor stop direction] = Terminate(tnew,coesout,JD0)
% This function specifies the event at which ode45 terminates.
muearth = 398600;

h = coesout(1);
ecc = coesout(2);
theta = coesout(3);
RAAN = coesout(4);
inc = coesout(5);
per = coesout(6);

[R,V] = coes2rvKAMME(h,inc,RAAN,ecc,per,theta,muearth);
r = norm(R);

alt = norm(r-6378);
lookfor = alt - 100; % = 0 when altitude = 100 km
stop = 1; % 1 means terminate at lookfor = 0; Otherwise 0
direction = -1; % -1 means zero crossing is from above
end 

function [xmax,imax,xmin,imin] = extrema(x)
%EXTREMA   Gets the global extrema points from a time series.
%   [XMAX,IMAX,XMIN,IMIN] = EXTREMA(X) returns the global minima and maxima 
%   points of the vector X ignoring NaN's, where
%    XMAX - maxima points in descending order
%    IMAX - indexes of the XMAX
%    XMIN - minima points in descending order
%    IMIN - indexes of the XMIN
%
%   DEFINITION (from http://en.wikipedia.org/wiki/Maxima_and_minima):
%   In mathematics, maxima and minima, also known as extrema, are points in
%   the domain of a function at which the function takes a largest value
%   (maximum) or smallest value (minimum), either within a given
%   neighbourhood (local extrema) or on the function domain in its entirety
%   (global extrema).
%
%   Example:
%      x = 2*pi*linspace(-1,1);
%      y = cos(x) - 0.5 + 0.5*rand(size(x)); y(40:45) = 1.85; y(50:53)=NaN;
%      [ymax,imax,ymin,imin] = extrema(y);
%      plot(x,y,x(imax),ymax,'g.',x(imin),ymin,'r.')
%
%   See also EXTREMA2, MAX, MIN
%   Written by
%   Lic. on Physics Carlos Adri?n Vargas Aguilera
%   Physical Oceanography MS candidate
%   UNIVERSIDAD DE GUADALAJARA 
%   Mexico, 2004
%
%   nubeobscura@hotmail.com
% From       : http://www.mathworks.com/matlabcentral/fileexchange
% File ID    : 12275
% Submited at: 2006-09-14
% 2006-11-11 : English translation from spanish. 
% 2006-11-17 : Accept NaN's.
% 2007-04-09 : Change name to MAXIMA, and definition added.
xmax = [];
imax = [];
xmin = [];
imin = [];
% Vector input?
Nt = numel(x);
if Nt ~= length(x)
 error('Entry must be a vector.')
end
% NaN's:
inan = find(isnan(x));
indx = 1:Nt;
if ~isempty(inan)
 indx(inan) = [];
 x(inan) = [];
 Nt = length(x);
end
% Difference between subsequent elements:
dx = diff(x);
% Is an horizontal line?
if ~any(dx)
 return
end
% Flat peaks? Put the middle element:
a = find(dx~=0);              % Indexes where x changes
lm = find(diff(a)~=1) + 1;    % Indexes where a do not changes
d = a(lm) - a(lm-1);          % Number of elements in the flat peak
a(lm) = a(lm) - floor(d/2);   % Save middle elements
a(end+1) = Nt;
% Peaks?
xa  = x(a);             % Serie without flat peaks
b = (diff(xa) > 0);     % 1  =>  positive slopes (minima begin)  
                        % 0  =>  negative slopes (maxima begin)
xb  = diff(b);          % -1 =>  maxima indexes (but one) 
                        % +1 =>  minima indexes (but one)
imax = find(xb == -1) + 1; % maxima indexes
imin = find(xb == +1) + 1; % minima indexes
imax = a(imax);
imin = a(imin);
nmaxi = length(imax);
nmini = length(imin);                
% Maximum or minumim on a flat peak at the ends?
if (nmaxi==0) && (nmini==0)
 if x(1) > x(Nt)
  xmax = x(1);
  imax = indx(1);
  xmin = x(Nt);
  imin = indx(Nt);
 elseif x(1) < x(Nt)
  xmax = x(Nt);
  imax = indx(Nt);
  xmin = x(1);
  imin = indx(1);
 end
 return
end
% Maximum or minumim at the ends?
if (nmaxi==0) 
 imax(1:2) = [1 Nt];
elseif (nmini==0)
 imin(1:2) = [1 Nt];
else
 if imax(1) < imin(1)
  imin(2:nmini+1) = imin;
  imin(1) = 1;
 else
  imax(2:nmaxi+1) = imax;
  imax(1) = 1;
 end
 if imax(end) > imin(end)
  imin(end+1) = Nt;
 else
  imax(end+1) = Nt;
 end
end
xmax = x(imax);
xmin = x(imin);
% NaN's:
if ~isempty(inan)
 imax = indx(imax);
 imin = indx(imin);
end
% Same size as x:
imax = reshape(imax,size(xmax));
imin = reshape(imin,size(xmin));
% Descending order:
[temp,inmax] = sort(-xmax); clear temp
xmax = xmax(inmax);
imax = imax(inmax);
[xmin,inmin] = sort(xmin);
imin = imin(inmin);
end