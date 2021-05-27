function [t0] = t0_solver(b,g,Cv,lcg,pw,u,DR,mb) 

% solve for SS trim and wetted waterline using savitsky methods
%Solver for lw (lp = lcg)
F = @(lw) lw*b*(.75-(5.21*Cv.^2/lw.^2+2.39).^(-1))-lcg; %lw - non dimensional wetted water length
x0=1.7;
[lw,fval]=fsolve(F,x0);

% Solver for Cl0 using Savitsky SS methods
Clb=mb*g/(.5*pw*u^2*b^2);
F = @(Cl0) Clb-Cl0+.0065*DR*Cl0.^.6;
x0=.14;
[Cl0,fval]=fsolve(F,x0);

% % SS estimate of tau
%lw=3;
F = @(t0) -Cl0+(t0).^1.1*(.012*lw.^.5+.0055*lw.^2.5/Cv^2);
x0=5;
%options=optimoptions('fsolve','Display','iter');
[t0,fval]=fsolve(F,x0);

%%
F = @(x) .89+(1.4/x)*(1-x^2)^.5;
x0=10;
%options=optimoptions('fsolve','Display','iter');
[x,fval]=fsolve(F,x0);

