%routine to establish jacobian, using interpolation on all forces, fast routine (MIMO)

function [dx,dx_c] = lin_interpsys4(e,XX,YY,tri,lcg,vcg,DR,ma,mb,A,B,t0,U0,int,I) 
g=32.2;
F=e(1);
S=e(2);
u=e(3);
v=e(4);
p=e(5);
r=e(6);
XI=[DR t0 v r F];

U_x=U0(1);
U_m=U0(2);
U_r=U0(3);

%% interpolation

[Yint]= interp_point(XX,YY,tri,XI,int); %XI=[DR t0 v r F]

%% Rudder modeling
z_rud=vcg; %ft

if I==1  %Actuator input: 1-mass, 2-rudder
X_rud=0;
Y_rud=0;
K_rud=0;
N_rud=0;
end

if I==2
%Simplified rudder model
SL=atand(v/u);
AoA=U_r-SL;
L=.5*(.0018)*u^2*AoA;
Y_rud=L;
X_rud=-.5*(4.6e-5*AoA^2+1.38e-4*abs(AoA)+9.2e-4)*u^2;
K_rud=-Y_rud*z_rud;
N_rud=-Y_rud*lcg;
end

C_X=U_x+X_rud; %setting thrust = drag for constant velocity
C_Y=Y_rud; %sway from rudder
C_K=ma*g*U_m+K_rud; % roll force from control mass & steering mechanism
% C_N=C_X*(U_m)*ma/(ma+mb)+N_rud;  %yaw force from steering helm + drag differential
C_N=N_rud; %thrust yaw coupling likely double counted in HD forces

F_hd=Yint'; 
%K=l_p*L_p-l_s*L_s; F_hd(3,:)=K; %Judge model incorporated
F_co=[C_X;C_Y;C_K;C_N];
loa=50/12; Iy=mb*(.25*loa)^2;
%  xg=0; q=0;yg=(U_m)*ma/(ma+mb);
 yg=0; xg=0; q=0;

%    u v   p         r
Cor=[0 0 0 -mb*v*(pi/180);        %X
    0 0 0 mb*u*(pi/180);        %Y
    0 0 0 0;         %K
    0 0 0 0];  %N

% RM - not needed since xdot and ydot are not part of the system  

% 6x6 system
e=[F;S;u;v;p;r];
dx=[zeros(2,4) eye(2); zeros(4,2) -A\B]*e+[zeros(2,6); zeros(4,2) -A\Cor]*e+[zeros(2,1);A\F_hd]; %full Cor
dx_c=[zeros(2,1);A\F_co];