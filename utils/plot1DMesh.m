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
% Description:  Plot 1D mesh with nodes and integration points						 	   %
%                                                                              %
% ---------------------------------------------------------------------------- %
function plot1DMesh()

  addpath('../src');

  in = Input;
	in.domain = [1,2];
	in.Ne     = 4;
	in.order  = 5;
	in.intrule = "GL";
	in.bctype  = {"periodic"};

	% Switch from periodic to dirichlet for plotting 
	% and dont plot left (right) part on the left (right) interval endopoint
	[el, fa] = CreateMesh1D(in);


	figure('Name',"Elements + quadrature + basis",'Position',[0,0,800,400]);
	hold on;

	Ne = el(1).Ne;

	for e = 1:Ne

		% Plot elements
		plot(el(e).nodes,[0,0],'-ko','MarkerSize',10,'MarkerFaceColor','k')

		% Plot quadrature points
		scatter(el(e).qp,zeros(1,length(el(e).qp)), 80,'o','filled')

		% The grid of quadrature points is too coarse
		% for i = 1:el(e).ndof
		% 	plot(el(e).qp,el(e).phi(i,:))
		% end
		
		x = linspace(el(e).nodes(1), el(e).nodes(2), 100);
		[P, dP] = LegendrePolynomial( x, in.order, el(e).nodes(1), el(e).nodes(2) );
		plot(x,P,'-k')

		ylim([1.1*min(min(P)), 1.1*max(max(P))])

	end

	xlim([in.domain(1) - 0.1, in.domain(2) + 0.1])


	figure('Name',"Faces and elements",'Position',[0,600,800,400]);
	hold on;

	Nf = fa(1).Nf;

	for	f = 1:Nf

		% Plot faces
		scatter(fa(f).node,0.1*f, 80,'o','filled')

		% Left element in blue
		plot(el( fa(f).elnr1 ).nodes,0.1*f*ones(1,el(e).nf), '-b');
		% Right element in red
		plot(el( fa(f).elnr2 ).nodes,0.1*f*ones(1,el(e).nf), '-r');

	end

	xlim([in.domain(1) - 0.1, in.domain(2) + 0.1])
	ylim([0,0.1*(Nf+1)])


	figure('Name',"Elements and faces",'Position',[900,0,800,400]);
	hold on;

	for e = 1:Ne

		% Plot elements
		plot([el(e).nodes(1),el(e).nodes(2)],0.1*e*ones(1,2),'-ko','MarkerSize',10,'MarkerFaceColor','k')

		% pLot centers
		scatter(el(e).center,0.1*e, 80,'o')

		% Left face in blue
		plot(fa( el(e).faces(1) ).node*ones(1,2), 0.1*e + 0.04*[-1,1],'-b')
		% Right face in red
		plot(fa( el(e).faces(2) ).node*ones(1,2), 0.1*e + 0.04*[-1,1],'-r')

		% Plot basis functions on elements
		% The grid of quadrature points is too coarse
		% for i = 1:el(e).ndof
		% 	plot(el(e).qp, 0.1*e + 0.04*el(e).phi(i,:))
		% end

		x = linspace(el(e).nodes(1), el(e).nodes(2), 50);
		P = LegendrePolynomial( x, in.order, el(e).nodes(1), el(e).nodes(2) );
		for i = 1:el(e).ndof
			plot(x,0.1*e + 0.04*P(i,:),'color', (i-1)*[1,1,1]./(el(e).ndof),'LineWidth',5)
		end

		% Plot basis functions on element faces
		phi = [ fa( el(e).faces(1) ).phi2'  ; ...
						fa( el(e).faces(2) ).phi1' ];

		x   = [ fa( el(e).faces(1) ).node  ; ...
					  fa( el(e).faces(2) ).node ];
		
		% Scale basis for plotting
		phi = 0.1*e + 0.04*phi; 

		for i = 1:el(e).ndof
			scatter(x,phi(:,i),(el(e).ndof+1-i)*200,'o','filled','MarkerFaceColor', (i-1)*[1,1,1]./(el(e).ndof))
		end

	end

	xlim([in.domain(1) - 0.1, in.domain(2) + 0.1])
	ylim([0,0.1*(Ne+1)])

end