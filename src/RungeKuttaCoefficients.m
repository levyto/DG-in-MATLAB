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
% Description:  Compute the coefficients of the Runge-Kutta method to be used  %
%               in the Solve function. The output is matrix A and vector b but %
%               for explicit SSP RK methods, these are not equal to the arrays %
%               in the Butcher's tableau. These are adjusted based on their    %
%               general formula in order to fit in the same framework.         %
%                                                                              %
% Input:   method ... Name of the method                                       %
%                                                                              %
% Output:       A ... A matrix of the RK method                                %
%               b ... b vector of the RK method                                %
%                                                                              %
% ---------------------------------------------------------------------------- %
function [A, b] = RungeKuttaCoefficients(method)

  if strcmp(method,"ForwardEuler")

    A = 1;
    b = 1;
    return;

  end

  if strcmp(method,"SSPRK2")
    A = [  1,   0;...
         1/2, 1/2];
    b = [  1; 1/2];
    return;
  end

  if strcmp(method,"SSPRK3")
    A = [  1,   0,   0;...
         3/4, 1/4,   0;...
         1/3,   0, 2/3];
    b = [  1; 1/4; 2/3];
    return;
  end

  if strcmp(method(1:11),"SSPRKLinear")

    % SSP RK methods for linear constant coefficients problems. Implemented
    % recursively as per: 
    %   [1] Gottlieb, S., Shu, C.-W., Tadmor, E.: Strong Stability-Preserving
    %       High-Order Time Discretization Methods, SIAM (2001)
    %
    % The second and third order method is equal to forward Euler and 2nd order
    % SSP RK method for nonlinear problems.
    %
    % For numerically obtained CFL_L2 numbers (i.e. our CFL/(2p+1)) see Table
    % 2.2 in
    %   [2] Cockburn, B, Shu, C.-W.: Runge–Kutta Discontinuous Galerkin Methods 
    %       for Convection-Dominated Problems, Journal of Scientific Computing
    %       (2001)
    % i.e.
    %         Table 2.2. The CFLL2 Numbers for Polynomials of Degree k  
    %                       and RK Methods of Order n
    %
    %     k    |   0      1      2      3      4      5      6      7      8
    %   ------------------------------------------------------------------------
    %   n = 1  | 1.000    *      *      *      *      *      *      *      *
    %   n = 2  | 1.000  0.333    *      *      *      *      *      *      * 
    %   n = 3  | 1.256  0.409  0.209  0.130  0.089  0.066  0.051  0.040  0.033
    %   n = 4  | 1.392  0.464  0.235  0.145  0.100  0.073  0.056  0.045  0.037
    %   n = 5  | 1.608  0.534  0.271  0.167  0.115  0.085  0.065  0.052  0.042
    %   n = 6  | 1.776  0.592  0.300  0.185  0.127  0.093  0.072  0.057  0.047
    %   n = 7  | 1.977  0.659  0.333  0.206  0.142  0.104  0.080  0.064  0.052
    %   n = 8  | 2.156  0.718  0.364  0.225  0.154  0.114  0.087  0.070  0.057
    %   n = 9  | 2.350  0.783  0.396  0.245  0.168  0.124  0.095  0.076  0.062
    %   n = 10 | 2.534  0.844  0.428  0.264  0.182  0.134  0.103  0.082  0.067
    %   n = 11 | 2.725  0.908  0.460  0.284  0.195  0.144  0.111  0.088  0.072
    %   n = 12 | 2.911  0.970  0.491  0.303  0.209  0.153  0.118  0.094  0.077 
    %
    % - Note that 1/(2p+1) is safe for polynomial of degree k and RK method of
    %   order n.
    % - The symbol '*' indicates that "the method is unstable when the ratio 
    %   dt/dx is held constant".

    m = str2num(method(12:end));

    if     m == 1
      [A,b] = RungeKuttaCoefficients("ForwardEuler");
      return;
    elseif m == 2
      [A,b] = RungeKuttaCoefficients("SSPRK2");
    else 
      A = diag(diag(ones(m,m)));
      b = ones(m,1);

      A(m,m) = 1/factorial(m);
      b(m)   = A(m,m);

      method_1 = method;
      method_1(12:end) = num2str(m-1);
      if (m-1) == 9
        method_1 = method_1(1:end-1);
      end
      A_1 = RungeKuttaCoefficients(method_1);
      for k = 1:(m-1)-1
        
        A(m,k+1) = A_1(m-1,k) / k;

      end

      A(m,1) = 1 - sum(A(m,2:m));

    end

  end

end