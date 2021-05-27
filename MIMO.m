% Linearizes system w/ 4th input for roll disturbance  5-13-21
% function [AA,BB,KK,K2] = SISO_LQR(x0,XX,YY,tri,lcg,vcg,DR,ma,mb,A,B,t0,U0,I,int) 
function [AA,BB] = MIMO(x0,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U0,I,int,CL)  % w 5DOF actuator dynamics 5-1
%baseline dynamics about setpoint
e=x0(3:10,1);
dD=.1;
% [dx,dx_c] = lin_interpsys3(e,XX,X,Y,K,N,Iy,lcg,vcg,DR,ma,mb,A,B,t0,U);
[dx,dx_c] = lin_interpsys8(e,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U0,int,I,CL);

%derivatives for system dynamics (AA matrix)
AA=zeros(8,8);
for j=1:8
e2=e; e2(j)=e(j)+dD;
[dx2,~] = lin_interpsys8(e2,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U0,int,I,CL);
AA(:,j)=(dx2-dx)/dD;
end

U2=U0; U2(2)=U0(2)+dD; % SI - mass
[~,dx_c2] = lin_interpsys8(e,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U2,int,I,CL);
BB(:,1)=(dx_c2-dx_c)./dD;

U2=U0; U2(3)=U0(3)+dD; % SI - rudder
[~,dx_c2] = lin_interpsys8(e,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U2,int,I,CL);
BB(:,2)=(dx_c2-dx_c)./dD;

U2=U0; U2(4)=U0(4)+dD; % SI - roll disturbance
[~,dx_c2] = lin_interpsys8(e,XX,YY,tri,lcg,vcg,zm,DR,ma,mb,A,B,t0,U2,int,I,CL);
BB(:,3)=(dx_c2-dx_c)./dD;
