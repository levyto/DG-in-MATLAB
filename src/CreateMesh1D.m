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
% Description:  Construct 1D equispaced FE mesh.                               %
%               NOTE: There are parts of the code commented out and labeled    %
%                     % NOT IMPLEMENTED. These relate to the other type of     %
%                     boundary conditions than periodic ones 
%                                                                              %
% Input:    input ... instance of class Input holding necessary setting for    %
%                     the start of the run                                     %
%                                                                              %
% Output:      el ... Array of Ne structs having the following properties:     %
%                     Ne .......... Number of elements                         %
%                          (global) Call with el(1).Ne)                        %
%                     DoF ......... Total number of Dofs                       %
%                          (global) Call with el(1).DoF)                       %
%                     NIP ......... Total number of integration points         %
%                          (global) Call with el(1).NIP)                       %
%                     elnr ........ ID of the element                          %
%                     h ... ....... Lenght of the element                      %
%                     nodes ....... Vector of endpoints                        %
%                     center ...... Center of the element                      %
%                     order ....... Polynomial order                           %
%                     ndof ........ Number of degrees of freedom               %
%                     offset ...... The pointer to the beginning of the        %
%                                   location in global coefficient vector      % 
%                                   which corresponds to given element         %
%                     nf .......... Number of faces adjacent to element        %
%                     faces ....... Array of face numbers creating the         %
%                                   element boundary                           %
%                     order_int ... Order which will be integrated exactly by  %
%                                   the quadrature                             %
%                     qp .......... Vector of quadrature points                %
%                     qw .......... Vector of quadrature points                %
%                     nip ......... Number of integration points               %
%                     phi ......... ndof-by-nip matrix of basis functions      %
%                                   evaluated at the quadrature points         %
%                     dphi ........ ndof-by-nip matrix of basis functions      %
%                                   derivatives evaluated at the quadr. pts    %
%                     mass ........ Mass matrix on element                     %
%                     inv_mass .... Mass matrix inverse on element             %
%                     sol ......... Solution reconstructed on the int. points  %
%                                                                              %
%              fa ... Array of Nf structs having the following properties:     %
%                     Nf .......... Number of faces                            %
%                          (global) Call with fa(1).Nf                         %
%                     facenr ...... ID of the face                             %
%                     node ........ x-coordinate of the face                   %
%                     phi1 ........ Basis function evaluated of the left       %
%                                   element evaluated at the face              %
%                     phi2 ........ Basis function evaluated of the right      %
%                                   element evaluated at the face              %
%                     n1 .......... Outward normal of element on the left      %
%                     n2 .......... Outward normal of element on the right     %
%                     elnr1 ....... Index of element on the left               %
%                     elnr2 ....... Index of element on the right              %
%                     bcnr ........ Type of boundary condition                 %
%                     sol ......... Trace of the solution on face (HDG)        %
%                     sol1 ........ Solution on the elem. on the on the left   %
%                     sol2 ........ Solution on the elem. on the on the right  %
%                                                                              %
% ---------------------------------------------------------------------------- %
function [el, fa] = CreateMesh1D(input)
  
  DoF = 0;
  NIP = 0;

  % ELEMENTS
  % --------
  for e = 1:input.Ne
    
    el(e).elnr   = e;
    el(e).h      = abs( input.domain(2) - input.domain(1) ) / input.Ne;
    el(e).nodes  = [ (e-1)*el(e).h, e*el(e).h ] + input.domain(1);
    el(e).center = mean(el(e).nodes);

    el(e).order  = input.order;
    el(e).ndof   = input.order + 1;
    el(e).offset = (e-1)*el(e).ndof+1;

    DoF = DoF + el(e).ndof;

    el(e).order_int = 2*input.order;  % integrand in mass matrix is of order 2p

    if     strcmp(input.intrule,"GL")
      [x, w] = GaussLegendreRule( el(e).order_int, el(e).nodes(1), el(e).nodes(2) );
    elseif strcmp(input.intrule,"LGL")
      [x, w] = LobattoGaussLegendreRule( el(e).order_int, el(e).nodes(1), el(e).nodes(2) );
    end
    el(e).qp  = x;
    el(e).qw  = w;
    el(e).nip = length(x);

    NIP = NIP + el(e).nip;

    % Compute basis functions
    [P, dP] = LegendrePolynomial( x, input.order, el(e).nodes(1), el(e).nodes(2) );
    el(e).phi  = P;
    el(e).dphi = dP;
    
    % Compute mass matrix
    el(e).mass     = ComputeMassMatrix( el(e).phi, el(e).qw );
    % el(e).mass = diag(diag(el(e).mass));
    el(e).inv_mass = inv(el(e).mass);

    el(e).sol = zeros(el(e).nip,1);

    % elements:    1     2
    %           |-----|-----|
    % faces:    1     2     3
    el(e).nf = 2;
    el(e).faces  = [ e, e+1 ];

    if e == input.Ne && strcmp(input.bctype,"periodic")
      el(e).faces(2) = 1;
    end

  end

  el(1).Ne = input.Ne;
  el(1).DoF = DoF;
  el(1).NIP = NIP;

  % Number of faces
% Nf = input.Ne+1;
  Nf = input.Ne; % ................................................. Periodicity
  fa(1).Nf = Nf;


  % Initialize faces
  fa( Nf ).facenr = Nf;

  for e = 1:input.Ne

    % FACES
    % -----
    % Visit each face only once - On each element only the left face but both on the last one
    visitF = 1;
%   if e == input.Ne                                                            % NOT IMPLEMENTED
%     visitF = el(e).nf;                                                        % NOT IMPLEMENTED
%   end                                                                         % NOT IMPLEMENTED
    for F = 1:visitF

      % Global face number
      f = el(e).faces(F); 
      fa(f).facenr = f;

      fa(f).node = el(e).nodes(F);

      fa(f).n1 = [ 1,0];
      fa(f).n2 = [-1,0];

      fa(f).elnr1 = f-1;
      fa(f).elnr2 = f;

      % Left boundary
      if f == 1 
        if     strcmp(input.bctype,"periodic")
          fa(f).elnr1 = input.Ne;
%       elseif strcmp(input.bctype,"Dirichlet")                                 % NOT IMPLEMENTED
%         fa(f).elnr1 = -1;                                                     % NOT IMPLEMENTED
        end
      end

%     % Right boundary                                                          % NOT IMPLEMENTED 
%     if f == Nf                                                                % NOT IMPLEMENTED  
%       if     strcmp(input.bctype, "periodic")                                 % NOT IMPLEMENTED 
%         fa(f).elnr2 = 1;                                                      % NOT IMPLEMENTED   
%       elseif strcmp(input.bctype, "Dirichlet")                                % NOT IMPLEMENTED   
%         fa(f).elnr2 = -1;                                                     % NOT IMPLEMENTED   
%       end                                                                     % NOT IMPLEMENTED   
%     end                                                                       % NOT IMPLEMENTED   

      if fa(f).elnr1 > 0
        P1 = LegendrePolynomial( el( fa(f).elnr1 ).nodes(2), ...
                                 input.order, ...
                                 el( fa(f).elnr1 ).nodes(1), ...
                                 el( fa(f).elnr1 ).nodes(2) );
        fa(f).phi1 = P1;
      else 
        fa(f).phi1 = NaN;
      end 

      if fa(f).elnr2 > 0
        P2 = LegendrePolynomial( el( fa(f).elnr2 ).nodes(1), ...
                                 input.order, ...
                                 el( fa(f).elnr2 ).nodes(1), ...
                                 el( fa(f).elnr2 ).nodes(2) );
        fa(f).phi2 = P2;
      else 
        fa(f).phi2 = NaN;
      end 

      fa(f).bcnr = 0;
      if     strcmp(input.bctype, "periodic")
        fa(f).bcnr = 0;
%     elseif strcmp(input.bctype, "Dirichlet")                                  % NOT IMPLEMENTED
%       fa(f).bcnr = 1;                                                         % NOT IMPLEMENTED
      end

      fa(f).sol = 0;

      fa(f).sol1 = 0;
      fa(f).sol2 = 0;

    end

  end

% if strcmp(bc,"periodic")                                                      % NOT IMPLEMENTED
%   % We need to keep the redundant face for plotting                           % NOT IMPLEMENTED
%   % and don't include it in assembly                                          % NOT IMPLEMENTED
%   fa(1).Nf = Nf-1;                                                            % NOT IMPLEMENTED
% end                                                                           % NOT IMPLEMENTED


end