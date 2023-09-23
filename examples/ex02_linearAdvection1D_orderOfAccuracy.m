% ---------------------------------------------------------------------------- %
% _______/\\\\\\\\\\\\_________________________/\\\\\\\\\\\\__________________ %
% _______\/\\\////////\\\_____________________/\\\//////////__________________ %
% ________\/\\\______\//\\\___________________/\\\____________________________ %
% _________\/\\\_______\/\\\__________________\/\\\____/\\\\\\\_______________ %
% __________\/\\\_______\/\\\ iscontinuous_____\/\\\___\/////\\\ alerkin______ %
% ___________\/\\\_______\/\\\__________________\/\\\_______\/\\\_____________ %
% ____________\/\\\_______/\\\___________________\/\\\_______\/\\\____________ %
% _____________\/\\\\\\\\\\\\/____________________\//\\\\\\\\\\\\/____________ %
% ______________\////////////_______________________\////////////_____________ %
%                                                                              %
% ---------------------------------------------------------------------------- %
%                                                                              %
% Description: Linear Scalar Advection Equation                                %
%                                                                              %
%              Solve                            u_t + (f(u))_x = 0, f(u) = au  %
%              with initial condition                   u(x,0) = u_0(x)        %
%              and periodic boundary conditions         u(0,t) = u(1,t)        %
%                                                                              %
%              Measure the order of accuracy of DK(p) method while using       %
%              Linear SSP Runge-Kutta method of order p+1.                     %
%                                                                              %
% ---------------------------------------------------------------------------- %
clear all; clc;
addpath('../src');
addpath('../models');

% Assign how many mesh refinements have to be done (> 1 to measure the order)
N_ELEMS = 2;
% Assign the highest polynomial approximation to be tested 
% The output will be results for p = O, 1, ..., N_ORDER
N_ORDER = 3;


function f = initCond(x)
  f = sin(2*pi*x);
end
in = Input;
in.model   = ScalarLinearAdvection1D( velocity = 1 );
in.domain  = [0,1];
in.Ne      = 16; % ............... Initial number of elements (will be doubled)
in.order   = 0; % ................ Initial order of polynomials
in.intrule = "GL";
in.timeint = "SSPRKLinear1"; % ... Intial RK method
in.CFL     = 0.95; % ............. To obtain reasonable results for p = 0 
in.endtime = 1; % ................ + domain + adv. velocity --> u(x,0) = u(x,1)
in.bctype  = {"periodic"};
in.ic      = @initCond;


% The error is computed as in
% [1] Lowrie, R. B.: Compact Higher-Order Numerical Methods For Hyperbolic
%     Conservation Laws, PhD. Thesis (1996)
% i.e.                          Ne
%      || w - w_h ||_2 = [ 1/Ne SUM ( P_j(w) - P_j(w_h) )^2 ]^(1/2)
%                               j=1
% where P_j(w), resp. P(w_h), is the mean of the analytic, resp. numerical, 
% solution at the left and right face of an element.
function err = L2NormError(el, fa, el0, fa0)
  Ne = el(1).Ne;
  err = 0;

  for e = 1:Ne

    sol = 0.5*( fa( el(e).faces(1) ).sol2 + fa( el(e).faces(2) ).sol1 );

    x1 = fa( el(e).faces(1) ).node;
    x2 = fa( el(e).faces(2) ).node;
    analytic = 0.5 * ( initCond(x1) + initCond(x2) );

    err = err + (sol - analytic)^2;

  end
  err = sqrt(err/Ne);

end


% Get the error for each Ne and p
L2Err = zeros(N_ORDER, N_ELEMS);
h     = zeros(1,N_ELEMS);
Ne0   = in.Ne

for j = 1:N_ORDER+1

  for k = 1:N_ELEMS

  % Initialize solver setup
    InitSolver(in);
  % Create 1D mesh
    [el, fa] = CreateMesh1D(in);
  % Initialize solution
    [coeffs, el, fa] = InitRun(el, fa, in);
  % Solve
    [coeffs, el, fa] = Solve(coeffs, el, fa, in);
  % Compute error
    L2Err(j,k) = L2NormError(el,fa);
    h(k) = 1/in.Ne;
  % Double the number of elements
    in.Ne = 2*in.Ne;

  end

% Reset number of elements
  in.Ne  = Ne0;
% Increase order
  in.order = in.order+1;
% Change time-integration method to be in.order + 1
  in.timeint(end) = num2str(j+1);

end 

% Plot errors in log-log plot
% ------------------------------------------------------------------------------
leg = "p = 0";
for j = 1:N_ORDER+1
  loglog(h,L2Err(j,:),'-o')
  hold on
  leg = [ leg; sprintf("p = %d",j)];
end
legend(leg);

% Print table
% ------------------------------------------------------------------------------
% Line
fprintf(" %s\n",repmat("-",1,78));

% Header
fprintf("  p\t\| ");
for j = 0:N_ORDER
  fprintf("%d\t",j);
end
fprintf("\n");

% Line
fprintf(" %s\n",repmat("-",1,78));

% First row
fprintf("  1/%d\t| ",1/h(1));
for j = 1:N_ORDER
  fprintf("-\t");
end 
fprintf("\n");

% Print measured orders
for k = 1:N_ELEMS-1
  fprintf("  1/%d\t| ",1/h(k+1));
  for j = 1:N_ORDER+1
    fprintf("%1.3f\t",log2( L2Err(j,k) / L2Err(j,k+1) ));
  end 
  fprintf("\n");
end

% Line
fprintf(" %s\n",repmat("-",1,78));
