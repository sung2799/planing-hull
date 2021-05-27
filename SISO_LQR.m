% function [AA,BB,KK,K2] = SISO_LQR(x0,XX,YY,tri,lcg,vcg,DR,ma,mb,A,B,t0,U0,I,int) 
function [AA,BB] = SISO_LQR(x0,XX,YY,tri,lcg,vcg,DR,ma,mb,A,B,t0,U0,I,int) 
%baseline dynamics about setpoint
e=x0(3:8,1);
dD=.1;
[dx,dx_c] = lin_interpsys4(e,XX,YY,tri,lcg,vcg,DR,ma,mb,A,B,t0,U0,int,I);

%Jacobian for system dynamics (AA matrix)
AA=zeros(6,6);
for j=1:6
e2=e; e2(j)=e(j)+dD;
[dx2,~] = lin_interpsys4(e2,XX,YY,tri,lcg,vcg,DR,ma,mb,A,B,t0,U0,int,I);
AA(:,j)=(dx2-dx)/dD;
end

%Jacobian for control dynamics (BB matrix)
if I==1
U2=U0; U2(2)=U0(2)+dD; % SI - mass
end
if I==2
U2=U0; U2(3)=U0(3)+dD; % SI - rudder
end
[~,dx_c2] = lin_interpsys4(e,XX,YY,tri,lcg,vcg,DR,ma,mb,A,B,t0,U2,int,I);
BB(:,1)=(dx_c2-dx_c)./dD;
