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
% Description:  Reconstruct the solution on every element and face by          %
%               evaluating the linear combination of basis functions           %
%                                                                              %
% Input:       el ... Array of elements to store the reconstructed solution    %
%              fa ... Array of    faces to store the reconstructed solution    %
%          coeffs ... Global coefficient vector                                %
%                                                                              %
% Output:      el ... Array of elements with stored reconstructed solution     %
%              fa ... Array of    faces with stored reconstructed solution     %
%                                                                              %
% ---------------------------------------------------------------------------- %
function [el, fa] = ReconstructSolution(el, fa, coeffs);

  Ne = el(1).Ne;

  for e = 1:Ne

    el(e).sol = zeros(size(el(e).qp))';

    for k = 1:el(e).nip 
      for i = 1:el(e).ndof
       el(e).sol(k) = el(e).sol(k) + coeffs( el(e).offset + i-1 ) * el(e).phi(i,k);
      end
    end

    % Use matrix-vector product
    % el(e).sol = ( el(e).phi )' * coeffs( el(e).offset : el(e).offset + el(e).ndof-1 );

  end

  Nf = fa(1).Nf;

  for f = 1:Nf

    fa(f).sol = 0; % TODO: HDG 

    elnr1 = fa(f).elnr1;
    if elnr1 > 0
      fa(f).sol1 = ( fa(f).phi1 )' * coeffs( el( elnr1 ).offset : el( elnr1 ).offset + el( elnr1 ).ndof-1 );
%   else                                                                        % NOT IMPLEMENTED
%     fa(f).sol1 = 0; % TODO: Dirichlet                                         % NOT IMPLEMENTED
    end

    elnr2 = fa(f).elnr2;
    if elnr2 > 0
      fa(f).sol2 = ( fa(f).phi2 )' * coeffs( el( elnr2 ).offset : el( elnr2 ).offset + el( elnr2 ).ndof-1 );
%   else                                                                        % NOT IMPLEMENTED
%     fa(f).sol2 = 0; % TODO: Dirichlet                                         % NOT IMPLEMENTED
    end

  end

end
