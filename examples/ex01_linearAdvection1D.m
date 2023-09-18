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
% ---------------------------------------------------------------------------- %
clear all; clc;
addpath('../src');
addpath('../models');


% Setup 1
% ------------------------------------------------------------------------------
in = Input;
in.model   = ScalarLinearAdvection1D( velocity = 2 ); % ... Advection velocity
in.domain  = [-1,1]; % .................................... Space domain
in.Ne      = 16; % ........................................ Number of elements
in.order   = 1; % ......................................... Polynomial order
in.intrule = "GL"; % ...................................... Quadrature rule
in.timeint = "ForwardEuler"; % ............................ Time integration
in.CFL     = 1; % ......................................... CFL number
in.endtime = 1; % ......................................... End time
in.bctype  = {"periodic"}; % .............................. Type of BCs
in.ic      = @(x) sin(pi*x); % ............................ Initial condition
in.showic  = 1; % ......................................... Plot IC and pause
in.showsol = 1; % ......................................... Plot solution on the fly


% Setup 2
% ------------------------------------------------------------------------------
% in = Input;
% in.model   = ScalarLinearAdvection1D( velocity = 2 );
% in.domain  = [-1,1];
% in.Ne      = 16;
% in.order   = 4;
% in.intrule = "GL";
% in.timeint = "SSPRK3";
% in.CFL     = 0.8;
% in.endtime = 1;
% in.bctype  = {"periodic"};
% in.ic      = @(x) exp(-64*x.^2);


% Setup 3
% ------------------------------------------------------------------------------
% function f = discontinuity(x)
%   f = ones(1,length(x));
%   f(x < -0.5) = 0;
%   f(x >  0.5) = 0;
% end
% in = Input;
% in.model   = ScalarLinearAdvection1D( velocity = 2 );
% in.domain  = [-1,1];
% in.Ne      = 50;
% in.order   = 2;
% in.intrule = "GL";
% in.timeint = "SSPRK3";
% in.CFL     = 0.8;
% in.endtime = 1;
% in.bctype  = {"periodic"};
% in.ic      = @discontinuity;


% Setup 4
% See [1] Krivodonova L.: Limiters for high-order discontinuous Galerkin
%         methods, JCP 226 (2007)
% ------------------------------------------------------------------------------
% function f = Krivodonova(x)
%   delta = 0.005;
%   alpha = 25;
%   b     = 0.75;
%   z     = 0.15;
%   beta  = log(4) / (36 * delta ^ 2);
%   ex    = @(x,y) exp(- beta .* (x - y) .^ 2);
%   F     = @(x,y) sqrt(max(1 - alpha ^ 2 * (x - y) .^ 2, 0));
%   f = 0.0 * (x >= 0 & x < 0.1) + ...
%       0.5 / 6 * (ex(x,z - delta) + ex(x,z + delta) + 4 .* ex(x,z)) .* (x >= 0.1 & x < 0.2) + ...
%       0.0 * (x >= 0.2 & x < 0.3) + ... 
%       0.5 * (x >= 0.3 & x < 0.4) + ...
%       0.0 * (x >= 0.4 & x < 0.5) + ...
%       (0.5 - abs(10 * (x - 0.55))) .* (x >= 0.5 & x < 0.6) + ...
%       0.0 * (x >= 0.6 & x < 0.7) + ...
%       0.5 / 6 * (F(x,b - delta) + F(x,b + delta) + 4 .* F(x,b)) .* (x >= 0.7 & x < 0.8) + ...
%       0.0 * (x >= 0.8 & x < 0.9) + ...
%       0.0 * (x >= 0.9 & x < 1); 
% end
% in = Input;
% in.model   = ScalarLinearAdvection1D( velocity = 1 );
% in.domain  = [0,1];
% in.Ne      = 100;
% in.order   = 1;
% in.intrule = "GL";
% in.timeint = "SSPRK3";
% in.CFL     = 1.0;
% in.endtime = 1;
% in.bctype  = {"periodic"};
% in.ic      = @Krivodonova;

% Initialize solver setup
InitSolver(in);

% Create 1D mesh
[el, fa] = CreateMesh1D(in);

% Initialize solution
[coeffs, el, fa] = InitRun(el, fa, in);

% Solve
[coeffs, el, fa] = Solve(coeffs, el, fa, in);

% Plot
PlotSolution(el, fa, coeffs);


