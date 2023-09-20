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
%              Plot eigenvalues of the L matrix (linear operator) in           %
%                                                                              %
%                 dU                                                           %
%                ---- = R(U) = L*U                                             %
%                 dt                                                           %
%                                                                              %
%              in order to assess the stability of the semi-discretized method %
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
in.order   = 5; % ......................................... Polynomial order
in.intrule = "GL"; % ...................................... Quadrature rule
in.bctype  = {"periodic"}; % .............................. Type of BCs
in.ic      = @(x) sin(2*pi*x); % .......................... Initial condition
in.test    = 1; % ......................................... Flag for testing


% Initialize solver setup
InitSolver(in);

% Create 1D mesh
[el, fa] = CreateMesh1D(in);

% Initialize solution
[coeffs, el, fa] = InitRun(el, fa, in);

% Assembly
[L, M, D, F] = AssembleLinearOperator(el, fa, in.model);

Lambda = eig(L);

figure();
plot(Lambda,'.','MarkerSize',15);
grid on;
xlabel('Re');
ylabel('Im');

fprintf("  Maximum of spatial operator eigenvalues: max(Re{lambda}) = "); 
fprintf("%1.3e\n\n",max(real(Lambda)));
