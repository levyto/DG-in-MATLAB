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
% Description: Scalar Inviscid Burgers' Equation                               %
%                                                                              %
%              Solve                            u_t + (f(u))_x = 0             %
%              with flux                                     f(u) = 0.5*u^2    %
%              initial condition                        u(x,0) = u_0(x)        %
%              and periodic boundary conditions         u(0,t) = u(1,t)        %
%                                                                              %
% ---------------------------------------------------------------------------- %
clear all; clc;
addpath('../src');
addpath('../models');

% Setup 1
% -----------------------------------------------------------------------------
in = Input;
in.model   = ScalarInviscidBurgers1D;
in.domain  = [0,1];
in.Ne      = 100;
in.order   = 1;
in.intrule = "GL";
in.timeint = "SSPRK3";
in.CFL     = 0.9;
in.endtime = 1;
in.bctype  = {"periodic"};
in.ic      = @(x) exp( -(10*(x-1/4)).^2 );


% Setup 2
% ----------------------------------------------------------------------------- 
% in = Input;
% in.model   = ScalarInviscidBurgers1D;
% in.domain  = [0,1];
% in.Ne      = 100;
% in.order   = 1;
% in.intrule = "GL";
% in.CFL     = 0.5;
% in.timeint = "SSPRK3";
% in.endtime = 1;
% in.bctype  = {"periodic"};
% in.ic      = @(x) 1/4 + 1/2*sin(pi*(2*x-1));

% Setup 3
% -----------------------------------------------------------------------------
% function f = discontinuity(x)
%   f = ones(1,length(x));
%   f(x < 0.1) = 0;
%   f(x > 0.3) = 0;
% end
% in = Input;
% in.model   = ScalarInviscidBurgers1D;
% in.domain  = [0,1];
% in.Ne      = 100;
% in.order   = 1;
% in.intrule = "GL";
% in.CFL     = 0.9;
% in.timeint = "SSPRK3";
% in.endtime = 1;
% in.bctype  = {"periodic"};
% in.ic      = @discontinuity;

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