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
% Description:  Compute L2 projection of the initial condition on the grid     %
%               points. When calling this function, ReconstructSolution must   %
%               follow.                                                        %
%                                                                              %
% Input:       el ... Array of elements                                        %
%              fa ... Array of faces                                           %
%          coeffs ... Global coefficient vector to be computed                 %
%              ic ... Function handle to evaluate the initial conditions       %
%                                                                              %
% Output:  coeffs ... Global coefficient vector to be computed                 %
%                                                                              %
% ---------------------------------------------------------------------------- %
function coeffs = InitialConditionL2Projection(el, fa, coeffs, ic);

  Ne = size(el,2);

  % Get solution coefficients on each element
  for e = 1:Ne

    % Asemble right-hand side
    rhs = zeros(el(e).ndof,1);
    for i = 1:el(e).ndof
      for k = 1:el(e).nip
        rhs(i) = rhs(i) + ic( el(e).qp(k) ) * el(e).phi(i,k) * el(e).qw(k);
      end
    end 

    % Solve linear system to obtain coefficients
    coeffs( el(e).offset : el(e).offset + el(e).ndof-1 ) = el(e).mass\rhs;

  end

end
