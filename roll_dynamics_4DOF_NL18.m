%5DOF: Surge(u),Sway(v),Roll(F),Yaw(S),Actuator, based on NL16 - 4/26 

function dy = roll_dynamics_4DOF_NL18(y,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U0,int,I,CL) % open loop

%Universal constants
pw=1.92; g=32.2; 
F=y(3);
S=y(4);
u=y(5);
v=y(6);
p=y(7);
r=y(8);
d1=y(9); %Actuator states (d1-position, d2-velocity)
d2=y(10);
XI=[DR t0 v r F];

U_r=U0(3); U_x=U0(1); U_p=U0(4);

if CL==0
    F_a=0;
end

if CL==1
F_a=U0(2); 
end

% actuator saturation
% if d1>.3
%     d1=.3;
% end
% if d1<-.3
%     d1=-.3;
% end

[Yint]= interp_point(XX,YY,tri,XI,int); %XI=[DR t0 v r F]

%% Rudder modeling
if I==1  %Actuator input: 1-mass, 2-rudder
X_rud=0;
Y_rud=0;
K_rud=0;
N_rud=0;
end

if I==2
z_rud=vcg; %ft

%Simplified rudder model
SL=atand(v/u);
AoA=U_r-SL;
L=.5*(.0018)*u^2*AoA;
Y_rud=L;
X_rud=-.5*(4.6e-5*AoA^2+1.38e-4*abs(AoA)+9.2e-4)*u^2;
K_rud=-Y_rud*z_rud;
N_rud=-Y_rud*lcg;

end

%% Forces
% C_X=U_x*cosd(U_r); %vector thrust 
% Cd=1e-5; A=(1/25)*loa*(.25); T=.5*Cd*pw*A*u^2; % drag dependent thrust

Ix=A(3,3)+ma*(d1^2+zm^2); %Dynamic roll moment of inertia
A2=A; A2(3,3)=Ix;

C_X=U_x+X_rud; %setting thrust = drag for constant velocity
C_Y=Y_rud-F_a; %sway from rudder
C_K=-F_a*zm+K_rud+U_p; % roll force from control mass & steering mechanism
% C_N=C_X*(U_m)*ma/(ma+mb)+N_rud;  %yaw force from steering helm + drag differential
C_N=N_rud; %thrust yaw coupling likely double counted in HD forces
C_d2=F_a;

F_hd=Yint'; 
%K=l_p*L_p-l_s*L_s; F_hd(3,:)=K; %Judge model incorporated
F_co=[C_X;C_Y;C_K;C_N];
loa=50/12; Iy=mb*(.25*loa)^2;
% yg=(U_m)*ma/(ma+mb)
yg=0; xg=0; q=0;

%Cor=[0; mb*r*(pi/180)*u; 0; 0];
%    u v   p         r
Cor=[0 0 mb*yg*q*(pi/180)^2 -mb*(xg*r*(pi/180)^2+v*(pi/180));        %X
    0 0 -mb*yg*p*(pi/180)^2 -mb*(yg*r*(pi/180)^2-u*(pi/180));        %Y
    -mb*yg*q*(pi/180) mb*yg*p*(pi/180) 0 -Iy*q*(pi/180)^2;         %K
    mb*(xg*r*(pi/180)+v) mb*(yg*r*(pi/180)-u) Iy*q*(pi/180)^2 0];  %N

RM=[cosd(S) -sind(S) 0 0;
    sind(S) cosd(S) 0 0;   
    0  0 1 0;           
    0 0 0 1];          

% Equation of motion
dy=[zeros(4) RM zeros(4,2); zeros(4) -A2\B zeros(4,2); zeros(2,10)]*y...     %Damping dynamics
    +[zeros(4,8) zeros(4,2);zeros(4) -A2\Cor zeros(4,2); zeros(2,10)]*y...   %Coriolis
    +[zeros(4,1);A2\(F_hd+F_co);0;ma\C_d2];                                  %External forces
% dy(7)=dy(7)+(Ix-ma*zm^2)^-1*(ma*g*d1*cosd(F)-ma*d1*p*(2*d2+zm*p));  %pd dynamics v1
dy(7)=dy(7)-ma/Ix*(2*d1*d2*p+zm*d1*p^2-g*d1*cosd(F));  %pd - actuator dynamics 
dy(9)=d2;   %controller velocity defined 
dy(10)=dy(10)+d1*p^2+g*sind(F)-zm*dy(7); %controller d dynamics
dy(5)=0; %set X direction acceleration to 0

if CL==0  %if Open Loop, freeze position of control mass
dy(10)=0;
end
warning('off')


