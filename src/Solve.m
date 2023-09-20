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
% Description:  Solution process driver to solve the problem on given mesh,    %
%               initial condition, and additional parameters given in input    %
%                                                                              %
% Input:       el ... Array of elements                                        %
%              fa ... Array of faces                                           %
%             sol ... Global coefficient vector to be solved for               %
%            nput ... Input parameters object                                  %
%                                                                              %
% Output:     sol ... Global coefficient vector of the solution                %
%              el ... Array of elements with updated solution                  %
%              fa ... Array of faces with updated solution                     %
%                                                                              %
% ---------------------------------------------------------------------------- %
function [sol, el, fa] = Solve(sol, el, fa, input)

  fprintf('  %s\n',repmat('·',1,76));
  fprintf('  #Timestep\tTime\t\tDeltaT\n');
  fprintf('  %s\n',repmat('·',1,76));
  tic

  Ne = el(1).Ne;
  
  dt = 1e10;
  t = 0;
  n = 0;
  while t < input.endtime

    % Compute time step size for RKDG
    for e = 1:Ne
      % Compute integral average of the solution
      solAvg = 0;
      for k = 1:el(e).nip
        solAvg = solAvg + el(e).qw(k) * el(e).sol(k);
      end
      solAvg = solAvg/el(e).h;
      % Take minimum over the elements
      dt = min(dt, input.CFL/(2*el(e).order+1) * el(e).h / input.model.EvalMaxEigenValue(solAvg) );
    end
    if t + dt > input.endtime - 1e-13
      dt = input.endtime - t;
    end

    n = n + 1;
    t = t + dt;
    fprintf("  %d\t\t%1.3e\t%1.3e\n",n,t,dt)

    % Runge-Kutta time integration
    [A, b] = RungeKuttaCoefficients(input.timeint);

    U = zeros(length(sol),length(b)+1);
    U(:,1) = sol;

    for s = 1:length(b)

      for q = 1:s
        U(:,s+1) = U(:,s+1) + A(s,q)*U(:,q);
      end

      U(:,s+1) = U(:,s+1) + b(s)*dt*Residual(el, fa, U(:,s), input.model);

      [el, fa] = ReconstructSolution(el, fa, U(:,s+1));

    end

    sol = U(:,end);

    % Plot solution during the run
    if input.showsol == 1

      hold off;
      PlotSolution(el, fa, sol);
      pause(1e-6);

    end

  end

  fprintf('  %s\n',repmat('·',1,76));
  fprintf("\n");
  fprintf("  Solved. ");
  toc
  fprintf("\n");


end