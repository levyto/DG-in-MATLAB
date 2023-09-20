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
% Description:  Plot the solution on each element 														 %
%                                                                              %
% Input:       el ... Array of elements                                        %
%              fa ... Array of faces                                           %
%          coeffs ... Global coefficient vector                                %
%                                                                              %
% ---------------------------------------------------------------------------- %
function PlotSolution(el, fa, coeffs, linecolor);

	Ne = size(el,2);

  for e = 1:Ne

  	if nargin == 3
  		linecolor = '-r';
  	end

	% Plot elements
	  % plot(el(e).nodes,[0,0],'-k')
	  % hold on;

	 % Plot solution on finer grid
		% x = linspace(el(e).nodes(1), el(e).nodes(2), 20);
	  % P = LegendrePolynomial( x, el(e).order, el(e).nodes(1), el(e).nodes(2) );
	  % sol = P' * coeffs( el(e).offset : el(e).offset + el(e).ndof-1 );
	  % plot(x,sol,'-k','LineWidth',2);

	 % Plot solution on quadrature points of the element
	  % plot(el(e).qp,el(e).sol,linecolor,'LineWidth',2)

	 % Plot solution on quadrature points and the end points of the element
	  if e == Ne % ................................................... Periodicity
	  	x = [ fa( el(e).faces(1) ).node, el(e).qp  , el(e).nodes(2) ];
	  	f = [ fa( el(e).faces(1) ).sol2, el(e).sol', fa( el(e).faces(2) ).sol1 ];
	  else
	  	x = [ fa( el(e).faces(1) ).node, el(e).qp  , fa( el(e).faces(2) ).node ];
	  	f = [ fa( el(e).faces(1) ).sol2, el(e).sol', fa( el(e).faces(2) ).sol1 ];
	  end

	  plot(x,f,linecolor,'LineWidth',2)
	  hold on;

	end

end