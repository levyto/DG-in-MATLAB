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
% Description:  Class for the user-defined input parameters. All parameters    %
%               are defined as NaNs in constructor and explicitly given by the %
%               user. The correctness of the inputs is performed in function   %
%               InitSolver().                                                  %
%                                                                              %
% ---------------------------------------------------------------------------- %
classdef Input

  properties

    % Mathematical model
    model % ...................... Object of model class stored in ../models 

    % Domain
    domain % ..................... Endpoints of the interval being discretized

    % Computational mesh 
    Ne % ......................... Number of elements

    % Polynomial approximation
    order % ...................... Order of the polynomials used on each element

    % Numerical integration
    intrule % .................... Quadrature rule for numerical integration

    % Time integration
    timeint % .................... Type of time-integration method
    timetype % ................... Explicit or implicit solution method
    CFL % ........................ CFL number
    endtime % .................... Final time of the simulation

    % Boundary conditions
    bctype % ..................... Array of types of boundary conditions
  % bc % ......................... Array of boundary condition function handles NOT IMPLEMENTED

    % Initial conditions
    ic % ......................... Initial condition function handle

    % Output
    showic % ..................... Plot initial condition and wait for input
    showsol % .................... Plot solution in every timestep

    % Test
    test % ....................... Flag for testing

  end
   
  methods

    function obj = Input()

      model   = NaN;
      domain  = NaN;
      Ne      = NaN;
      order   = NaN;
      intrule = NaN;
      timeint = NaN;
      timetype= NaN;
      CFL     = NaN;
      endTime = NaN;
      bctype  = {NaN};
      bc      = {@(x)NaN};
      ic      = @(x)NaN;
      showic  = NaN;
      showsol = NaN;
       
    end

  end

end