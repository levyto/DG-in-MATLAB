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
% Description:  Plot L2 projection of given function on given mesh defined on  %
%               [-1,1] and given polynomial space on every element             %
%                                                                              %
% ---------------------------------------------------------------------------- %
function plotL2Projection()

  addpath('../src');

  in = Input;
  in.domain  = [-1,1];
  in.order   = 2;
  in.Ne      = 5;
  in.intrule = "GL";
  in.bctype  = {"periodic"};

  k = 1;
  for i = 1:2
    for j = 1:2

      figure(k,'Position',[900*(i-1),600*(j-1),800,400]);
      hold on;

      L2ProjectionOfFun(i+j,in);
      
      k = k+1;

    end
  end

end


% ---------------------------------------------------------------------------- %
%                                                                              %
% Description:  Evaluate function on points x based on type                    %
%                                                                              %
% Input:        x ... Points where the function is to be evaluated             %
%              id ... Type of the function                                     %
%                                                                              %
% Output:       f ... Function evaluated at points x                           %
%            name ... Function expression                                      %
%                                                                              %
% ---------------------------------------------------------------------------- %
function [f, name] = fun(x,id)

  switch id
    case 1
      name = "sin(\pi x)";
      f = sin(pi*x);
    case 2
      name = "x^{10}";
      f = x.^10;
    case 3
      name = "e^{-64x^2}";
      f = exp(-64*x.^2);
    case 4
      name = "Discontinuous";
      f = zeros(1,length(x));
      f = f + 1;
      f(x > 0) = 0.125;
  end

end


% ---------------------------------------------------------------------------- %
%                                                                              %
% Description:  Compute L2 projection of given function on specified mesh.     %
%               Plot the result                                                %
%                                                                              %
% Input:        F ... Index of the function                                    %
%           input ... Input parameters object                                  %
%                                                                              %
% ---------------------------------------------------------------------------- %
function L2ProjectionOfFun(F, input)

  % F ...Function to be projected: 
  %   1 = sin 
  %   2 = 10th deg. poly.
  %   3 = Gaussian
  %   4 = Discontinuity inside element

  [el, fa] = CreateMesh1D(input);

  Ne  = el(1).Ne;
  DoF = el(1).DoF;

  % Initialize global coefficient vector
  coeffs = zeros(DoF,1);


  % Get solution coefficients on each element
  for e = 1:Ne

    % Asemble right-hand side
    rhs = zeros(el(e).ndof,1);
    for i = 1:el(e).ndof
      for k = 1:el(e).nip
        rhs(i) = rhs(i) + fun( el(e).qp(k), F ) * el(e).phi(i,k) * el(e).qw(k);
      end
    end 

    % Solve linear system to obtain coefficients
    coeffs( el(e).offset : el(e).offset + el(e).ndof-1 ) = el(e).mass\rhs;

  end

  % Reconstruct solution on every element
  [el, fa] = ReconstructSolution(el, fa, coeffs);

  Min = 1e3; Max = -1e3;
  for e = 1:Ne

    % Plot elements
    plot(el(e).nodes,[0,0],'-k')

    % Plot quadrature points
    scatter(el(e).qp,zeros(1,length(el(e).qp)), 80,'o','filled','MarkerFaceColor','b')

    % Plot the original function
    x = linspace(el(e).nodes(1), el(e).nodes(2), 100);
    f = fun(x,F);
    plot(x,f,'color', [211,211,211]./256,'LineWidth',10);

    % Plot reconstructed (and projected) function on finer grid
    P = LegendrePolynomial( x, input.order, el(e).nodes(1), el(e).nodes(2) );
    sol = P' * coeffs( el(e).offset : el(e).offset + el(e).ndof-1 );
    plot(x,sol,'color', [169,169,169]./256,'LineWidth',5);

    % Plot reconstructed (and projected) function on quadrature points
    plot(el(e).qp,el(e).sol,'-ro','MarkerSize',10,'MarkerFaceColor','r','LineWidth',2)

    % Plot reconstructed (and projected) function on faces
    Fsol = [ fa( el(e).faces(1) ).sol2, fa( el(e).faces(2) ).sol1 ];
    plot(el(e).nodes,Fsol,'ko','MarkerSize',10,'MarkerFaceColor',[169,169,169]./256,'LineWidth',2)


    % Plot element limits
    plot(el(e).nodes(1)*ones(1,2),[-1e3,1e3],'--k')
    if e == Ne
      plot(el(e).nodes(2)*ones(1,2),[-1e3,1e3],'--k')
    end

    Min = min( Min, min( min(f), min(sol) ) );
    Max = max( Max, max( max(f), max(sol) ) );

  end

  legend('elements',...
         'quadrature points x_q^{(e)}',...
         'original function',...
         'projected function',...
         'projected function on x_q^{(e)}',...
         'projected function on faces');

  [~,name] = fun(0,F);
  title(name)

  xlim([input.domain(1) - 0.1, input.domain(2) + 0.1])
  ylim([ min( -0.1, 1.1*Min ), 1.1*Max ])

end