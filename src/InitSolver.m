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
% Description:  Initialize the solver by checking the user-defined solver      %
%               parameters                                                     %
%                                                                              %
% Input:    input ... Input parameters object                                  %
%                                                                              %
% ---------------------------------------------------------------------------- %
function InitSolver(input)

  fprintf("\n");
  fprintf(" _______/\\\\\\\\\\\\\\\\\\\\\\\\_________________________/\\\\\\\\\\\\\\\\\\\\\\\\__________________ \n");
  fprintf(" _______\\/\\\\\\////////\\\\\\_____________________/\\\\\\//////////__________________ \n");
  fprintf(" ________\\/\\\\\\______\\//\\\\\\___________________/\\\\\\____________________________ \n");
  fprintf(" _________\\/\\\\\\_______\\/\\\\\\__________________\\/\\\\\\____/\\\\\\\\\\\\\\_______________ \n");
  fprintf(" __________\\/\\\\\\_______\\/\\\\\\ iscontinuous_____\\/\\\\\\___\\/////\\\\\\ alerkin______ \n");
  fprintf(" ___________\\/\\\\\\_______\\/\\\\\\__________________\\/\\\\\\_______\\/\\\\\\_____________ \n");
  fprintf(" ____________\\/\\\\\\_______/\\\\\\___________________\\/\\\\\\_______\\/\\\\\\____________ \n");
  fprintf(" _____________\\/\\\\\\\\\\\\\\\\\\\\\\\\/____________________\\//\\\\\\\\\\\\\\\\\\\\\\\\/____________ \n");
  fprintf(" ______________\\////////////_______________________\\////////////_____________ \n");
  fprintf("\n");
  
  if isequal(class(input.model),'ScalarLinearAdvection1D') || ...
     isequal(class(input.model),'ScalarInviscidBurgers1D')
    fprintf("  Solving:   %s\n",input.model.Name);
  else 
    msg = sprintf([...
          "Wrong 'Input.model'. Known types are:\n"...
          "\t- 'ScalarLinearAdvection1D'\n"...
          ]);
    error(msg);
  end
  fprintf("\n");

  % Testing
  if input.test == 1
    return;
  end


  % -----------------------------------
  % 1D Scalar Linear Advection Equation
  % 1D Scalar Invicid Burgers Equation
  % -----------------------------------
  if isequal(class(input.model),'ScalarLinearAdvection1D') || ...
     isequal(class(input.model),'ScalarInviscidBurgers1D')

    if isnan(input.domain)
      error("Domain not specified");
    elseif input.domain(2) <= input.domain(1)
      error("Incorrectly defined domain: [%.3f, %.3f]",input.domain(1),input.domain(2));
    else
      fprintf("%-25s[%.3f,%.3f]\n","  Computational domain:",input.domain(1),input.domain(2));
    end

    if isnan(input.Ne)
      error("Number of elements not specified");
    elseif input.Ne < 1
      error("Incorrectly specified number of elements: Input.Ne = %d",input.Ne);
    else
      fprintf("%-25s%d\n","  Number of elements:",input.Ne);
    end

    if isnan(input.order)
      error("Order of polynomials on each element not specified");
    elseif input.order < 0
      error("Incorrectly specified polynomial order: Input.order = %d",input.order);
    else
      fprintf("%-25s%d\n","  Order of polynomials:",input.order);
    end

    if isnan(input.intrule)
      error("Numerical quadrature rule not specified");
    elseif !(strcmp(input.intrule,"GL" ) || ...
             strcmp(input.intrule,"LGL"))
      msg = sprintf([...
            "Quadrature type not known. Known types are:\n"...
            "\t- 'GL'  \t(Gauss-Legendre)\n"...
            "\t- 'LGL' \t(Lobatto-Gauss-Legendre)\n"...
            ]);
      error(msg);
    else
      if     strcmp(input.intrule,"GL")
        rule = "Gauss-Legendre";
      elseif strcmp(input.intrule,"LGL")
        rule = "Lobatto-Gauss-Legendre";
      end
      fprintf("%-25s%s\n","  Integration rule:",rule);
    end

    if isnan(input.timeint)
      error("Time-integration method not specified");
    elseif !(strcmp(input.timeint,"ForwardEuler") || ...
             (length(input.timeint) > 5 && ...
              strcmp(input.timeint(1:5),"SSPRK"))       ) 
      msg = sprintf([...
            "Time-integration method not known. Known types are:\n"...
            "\t- 'ForwardEuler'\t(Explicit First-Order Forward Euler)\n"...
            "\t- 'SSPRK2'      \t(Explicit Two-Stage Second-Order SSP Runge-Kutta)\n"...
            "\t- 'SSPRK3'      \t(Explicit Two-Stage Third-Order SSP Runge-Kutta)\n"...
            "\t- 'SSPRKLinearX'\t(Explicit X-Stage X-Order Linear SSP Runge-Kutta)\n"...
            ]);
      error(msg);
    else
      timetype = "explicit";
      if     strcmp(input.timeint,"ForwardEuler")
        timeint = "Explicit forward Euler, 1st order";
      elseif strcmp(input.timeint,"SSPRK2")
        timeint = "Explicit SPP Runge-Kutta, 2nd order, 2 stages";
      elseif strcmp(input.timeint,"SSPRK3")
        timeint = "Explicit SPP Runge-Kutta, 3rd order, 3 stages";
      elseif length(input.timeint) > 11
        if !strcmp(input.timeint(1:11),"SSPRKLinear")
          error(msg);
        elseif str2num(input.timeint(12:end)) < 1
          error(msg);
        elseif str2num(input.timeint(12:end)) == 1
          timeint = "Explicit forward Euler, 1st order";
        elseif str2num(input.timeint(12:end)) == 2
          timeint = "Explicit SPP Runge-Kutta, 2nd order, 2 stages";
        else
          RKorder = str2num(input.timeint(12:end));
          timeint = ["Explicit SPP Linear Runge-Kutta, ",num2str(RKorder)];
          if RKorder == 3
            timeint = [timeint,"rd"];
          else
            timeint = [timeint,"th"];
          end
          timeint = [timeint," order, ",num2str(RKorder)," stages"];
        end
        fprintf("%-25s%s\n","  Time integration:",timeint);
    end


    if isnan(input.CFL)
      error("CFL number not specified");
    elseif input.CFL <= 0
      error("Incorrectly specified CFL number: Input.CFL = %d",input.CFL);
    else
      fprintf("%-25s%d\n","  CFL:",input.CFL);
    end

    if isnan(input.endtime)
      error("End time not specified");
    elseif input.endtime <= 0
      error("Incorrectly specified end time: Input.endtime = %d",input.endtime);
    else
      fprintf("%-25s%d\n","  End time:",input.endtime);
    end

    if isnan(input.bctype{1})
      error("Type of boundary conditions not specified");
    elseif !iscell(input.bctype)
      msg = sprintf([...
            "Incorrectly specified boundary conditions. "...
            "For e.g. periodic BCs use:\n"...
            "\t 'Input.bctype  = {\"periodic\"};'\n"
            ]);
      error(msg);
    elseif strcmp(input.bctype{1},"periodic") && ...
           size(input.bctype,2) > 1
      msg = sprintf([...
            "Incorrectly specified boundary conditions. "...
            "For periodic BCs use:\n"...
            "\t 'Input.bctype  = {\"periodic\"};'\n"
            ]);
      error(msg);
    elseif !(strcmp(input.bctype{1},"periodic" ) ) % || ...
%            strcmp(input.bctype{1},"Dirichlet") || ...                         % NOT IMPLEMENTED
%            strcmp(input.bctype{2},"Dirichlet"))                               % NOT IMPLEMENTED
      msg = sprintf([...
            "BC type not known. Known types are:\n"...
            "\t- 'periodic'\n"...
%           "\t- 'Dirichlet'\n"...                                              % NOT IMPLEMENTED
            ]);
      error(msg);
    else
      fprintf("%-25s","  Boundary conditions:");
      for i = 1:size(input.bctype,2)
        fprintf("%s",input.bctype{i});
        if i != size(input.bctype,2)
          fprintf(", ");
        end
      end
      fprintf("\n");
    end

    if isnan(input.ic(0))
      error("Initial conditions not specified");
    elseif !isa(input.ic,'function_handle')
      msg = sprintf([...
            "Incorrectly specified initial condition. "...
            "Specify IC as a function handle:\n"...
            "\t 'Input.ic  = @initialCondition;'\n"...
            "where initialCondition is user-defined function.\n"...
            "For simple IC like sin wave use:\n"...
            "\t 'Input.ic = @(x) sin(pi*x)'\n"
            ]);
      error(msg);
    end

  end
  % End 1D Scalar Linear Advection Equation
  % End 1D Scalar Invicid Burgers Equation

end