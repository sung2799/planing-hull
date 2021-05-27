% Req: XI, lcg, DR, sim_setup, & interp_model(lcg) for XX,YY,tri

function [Yint]= interp_point(XX,YY,tri,XI,int) %XI=[DR t0 v r F]

%%
% format for XI=[DR t0 v r F];

if int==1  % allows for use of MATLAB internal interpolator - much slower method
load('D_m.mat','XX_m','X_m','Y_m','K_m','N_m')
X = griddatan(XX_m,X_m,XI);
% Y=0; %Y model with linear damping matrix
Y = griddatan(XX_m,Y_m,XI); %Y test data yields unintuitive results
K = griddatan(XX_m,K_m,XI);
N = griddatan(XX_m,N_m,XI);
Yint=[X Y K N];
end

if int==2  % our custom script that only requires a single simplex calculation - much faster
[t,p] = tsearchn(XX,tri,XI);
YY_X=YY(:,1); YY_Y=YY(:,2); YY_K=YY(:,3); YY_N=YY(:,4);
X = p*YY_X(tri(t,:));
Y = p*YY_Y(tri(t,:));
K = p*YY_K(tri(t,:));
N = p*YY_N(tri(t,:));
Yint=[X Y K N];

end
