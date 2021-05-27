% Vessel particulars
pw=1.92; g=32.2;
mb=11.49/g;
loa=50/12;
b=9/12;
% vcg=5.29/12; %based on Judge
vcg=.25; %based B&K rudder data (K/Y)
% vcg=0; %Reduce vcg to reduce roll to fit Kim 2017 tests (3/23)
Cv=3;
u=Cv*(g*b)^.5;
lk=1.95*b;
[t0] = t0_solver(b,g,Cv,lcg,pw,u,DR,mb);
lc=lk-(b/pi)*tand(DR)/tand(t0);
lw=(lk+lc)/(2*b);
Cl=.7;
ma=.1*mb; %actuation mass used in active controller

% Mass matrices
U_m=d0;
Ix=mb*((.4*b)^2+(U_m*ma/(ma+mb))^2); %Roll: Using masumi 2019, similar to rectangular prism
% Iz=mb*((.25*loa)^2+(loa/2-lcg)^2);  %Yaw: parallel axis thm for shift in lcg
Iz=mb*((.2*loa)^2+(loa/2-lcg)^2);  %Reduced to better match Ad/TD data. 3-29
kB=.06641+.00716*DR+.0003861*DR^2; %hull added mass fxn
Yvd=-b^2*pw*tand(DR)*kB*(lk+2*lc)/12;
Nvd=-b^2*pw*tand(DR)*kB*(lk^2+2*lk*lc+3*lc^2)/48;
Nrd=-b^2*pw*tand(DR)*kB*(lk^3+2*lk^2*lc+3*lc^2*lk+4*lc^3)/120;
Kpd=-.010237*pw*b^5*lw*(1-sind(DR));
AR=[mb 0 0 0; 0 mb 0 0; 0 0 Ix -Ix/5;0 0 -Ix/5 Iz]; %Rigid body mass matrix (slugs-ft^2)
AM=-[.08*mb 0 0 0;  %added mass is (-) of maneuvering derivatives for acceleration
    0 Yvd 0 Nvd;
    0 0 Kpd 0;
    0 Nvd 0 Nrd];

A=AR+AM; %Mass(Rigid) + Mass(Added)

% Damping matrix
Kp=-(1-sind(DR))*(.029*Cv+.02*lw)*pw*g*b^4*(b/g)^.5;
% Kp=.2*Kp; %Reducing roll damping to achieve overshoot as shown in Kim 2017 (3/23)
B=-[0 0 0 0; 
    0 0 0 0; 
    0 0 Kp 0; 
    0 0 0 0]; %+damping is (-) of maneuvering derivatives for velocity
