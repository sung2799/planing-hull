%interpolated system, custom interpolator 1-18-21

function [XX,YY,tri]= interp_model(lcg)

%% Setup interpolation model
% D_b = xlsread('data_interp2.xls',2,'A2:J97'); %sheet 2 for body axes data
D_w = xlsread('data_interp2.xls',1,'A2:J177');  

% Establish measurement axes data

Xw=D_w(:,7); 
Yw=D_w(:,8);
Kw=D_w(:,9);
Nw=D_w(:,10);
SS=D_w(:,3);  %side slip angle

for j=1:length(D_w)
XYm(:,j)=[cosd(SS(j)) sind(SS(j)); -sind(SS(j)) cosd(SS(j))]\[Xw(j); Yw(j)];  %rotation matrix
Km(j)=Kw(j)/cosd(SS(j));
Nm(j)=Nw(j)-(-1.875+lcg)*Yw(j);
end
Xm=-XYm(1,:)'; Ym=XYm(2,:)'; Km=Km'; Nm=Nm'; %Note: Xw sign is backward
D_m=[D_w(:,1:6) Xm Ym Km Nm];

E_m=D_m; %Extrapolation Model (based on symmetric v,r,F,Y,K,N)
E_m(:,4)=-E_m(:,4); %v
E_m(:,5)=-E_m(:,5); %r
E_m(:,6)=-E_m(:,6); %F
E_m(:,8)=-E_m(:,8); %Y
E_m(:,9)=-E_m(:,9); %K
E_m(:,10)=-E_m(:,10); %N
D_m=[D_m;E_m];
XX=[D_m(:,1:2) D_m(:,4:6)]; %side slip angle removed since not needed
X=D_m(:,7); Y=D_m(:,8); K=D_m(:,9); N=D_m(:,10); 
XX_m=XX; K_m=K; Y_m=Y; X_m=X; N_m=N; %save for use in built-in call function
save('D_m.mat','XX_m','K_m','Y_m','X_m','N_m');
% YY=[X Y K N];  %Raw Force vector

[XX,~,ind2] = unique(XX,'rows');
YY_X = accumarray(ind2,X,[size(XX,1),1],@mean); % YY_X=X(ind1);
YY_Y = accumarray(ind2,Y,[size(XX,1),1],@mean); % YY_Y=Y(ind1);
YY_K = accumarray(ind2,K,[size(XX,1),1],@mean); % YY_K=K(ind1);
YY_N = accumarray(ind2,N,[size(XX,1),1],@mean); % YY_N=N(ind1);
YY=[YY_X,YY_Y,YY_K,YY_N];  %averaged force vector

tri = delaunayn(XX);
% % output 4 Y matrices for each individual force


