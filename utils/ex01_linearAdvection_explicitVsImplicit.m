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
%              and periodic boundary conditions         u(a,t) = u(b,t)        %
%                                                                              %
%              Compare explicit (Forward Euler) and implicit (Backward Euler   %
%              or Explicit firstage singly diagonally implicit Runge-Kutta     %
%              method of 2nd order) solution approach for the linear problem   %
%              given by:                                                       %
%                                                                              %
%                 dU                                                           %
%                ---- = R(U) = L*U                                             %
%                 dt                                                           %
%                                                                              %
% ---------------------------------------------------------------------------- %
clear all; clc;
addpath('../src');
addpath('../models');

% ==============================================================================
% Setup
% ==============================================================================
in = Input;
in.model   = ScalarLinearAdvection1D( velocity = 1 ); % ... Advection velocity
in.domain  = [0,1]; % ..................................... Space domain
in.Ne      = 50; % ........................................ Number of elements
in.order   = 1; % ......................................... Polynomial order
in.intrule = "GL"; % ...................................... Quadrature rule
in.timeint = "ForwardEuler"; % ............................ Time integration
in.CFL     = 0.6; % ....................................... CFL number
in.endtime = 1; % ......................................... End time
in.bctype  = {"periodic"}; % .............................. Type of BCs
in.ic      = @(x) sin(2*pi*x); % .......................... Initial condition
in.showic  = 0; % ......................................... Plot IC and pause
in.showsol = 1; % ......................................... Plot solution on the fly

implicit = "BackwardEuler"; % ........... 
% implicit = "ESDIRK2"; % ............... Type of implicit scheme
dt = 0.01; % ............................ Timestep for implicit time-integration
% ==============================================================================

% Initialize solver setup
InitSolver(in);

% Create 1D mesh
[el, fa] = CreateMesh1D(in);

% ------------------------------------------------------------------------------
% EXPLICIT
% ------------------------------------------------------------------------------
% tic
% % Initialize solution
% [coeffs, el, fa] = InitRun(el, fa, in);
% % Solve
% [coeffs, el, fa] = Solve(coeffs, el, fa, in);
% % Plot
% hold off;
% PlotSolution(el, fa, coeffs,'-b');

% pause(0.1);

% ------------------------------------------------------------------------------
% IMPLICIT
% ------------------------------------------------------------------------------
tic
% Initialize solution
[coeffs, el, fa] = InitRun(el, fa, in);

DoF = el(1).DoF;
t = 0;
n = 0;
while t < in.endtime

  if t + dt > in.endtime - 1e-13
    dt = in.endtime - t;
  end

  n = n + 1;
  t = t + dt;
  fprintf("  %d\t\t%1.3e\t%1.3e\n",n,t,dt);

  % Assembly
  I = eye(DoF);
  L = AssembleLinearOperator(el, fa, in.model);


  % Backward Euler -------------------------------------------------------------
  if strcmp(implicit,"BackwardEuler")

    % Solve system
    coeffs = ( I - dt*L ) \ coeffs;

    % Reconstruct solution on elements and faces
    [el, fa] = ReconstructSolution(el, fa, coeffs);


  % 2nd order three-stage ESDIRK -----------------------------------------------
  elseif strcmp(implicit,"ESDIRK2")

    % Butcher's A matrix
    A = [             0,             0,           0 ;...
            1-1/sqrt(2),   1-1/sqrt(2),           0 ;...
          1/(2*sqrt(2)), 1/(2*sqrt(2)), 1-1/sqrt(2) ]; 
    stages = size(A,1);

    U = zeros(DoF,stages);

    % First stage is explicit
    U(:,1) = coeffs;

    for s = 2:stages

      % Assemble RHS
      R = coeffs;
      for q = 1:s-1
        R = R + dt*A(s,q)*L*U(:,q);
      end 

      % Solve system
      U(:,s) = (I - dt*A(s,s)*L) \ R;

      % Reconstruct solution on elements and faces
      [el, fa] = ReconstructSolution(el, fa, U(:,s));

    end

    coeffs = U(:,end);

  end


  % Plot solution during the run
  if in.showsol == 1

    hold off;
    PlotSolution(el, fa, coeffs);
    pause(1e-6);
    % pause

  end

end

fprintf('  %s\n',repmat('·',1,76));
fprintf("\n");
fprintf("  Solved. ");
toc
fprintf("\n");


% Plot
PlotSolution(el, fa, coeffs,'r');


plt(1) = plot(NaN,NaN,'-b','LineWidth',2);
plt(2) = plot(NaN,NaN,'-r','LineWidth',2);
legend(plt, {'explicit', 'implicit'})
