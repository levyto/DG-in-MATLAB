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
%              Plot sparsity patterns of matrices appearing in the L matrix    %
%              (linear operator) in                                            %
%                                                                              %
%                 dU                                                           %
%                ---- = R(U) = L*U                                             %
%                 dt                                                           %
%                                                                              %
%              where L =  - M^{-1} ( -D + F )                                  %
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
in.Ne      = 5; % ......................................... Number of elements
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

figure()

subplot(2,2,1);
spy(M,'b.',15);
title('Mass matrix sparsity pattern');

subplot(2,2,2);
spy(D,'r.',15);
title('Differentiation matrix sparsity pattern');

subplot(2,2,3);
spy(F,'m.',15);
title('Flux matrix sparsity pattern');

subplot(2,2,4);
spy(L,'k.',15);
title('Global matrix sparsity pattern');