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
% Description:  Compute the reisudal vector of the fully-discretized problem   %
% Input:       el ... Array of elements                                        %
%              fa ... Array of faces                                           %
%             sol ... Global coefficient vector to be computed                 %
%           model ... Mathematical model defining the flux                     %
%                                                                              %
% Output:       R ... DG residual for given mathematical model                 %
%                                                                              %
% ---------------------------------------------------------------------------- %
function R = Residual(el, fa, sol, model)

  Ne = el(1).Ne;
  Nf = fa(1).Nf;

  % Residual vector
  R = zeros(size(sol));

  % Mass matrix
  M = zeros(length(sol));

  D = zeros(length(sol));


  % Ovel all elements
  for e = 1:Ne

    % Indices corresponding to element e in the global matrix
    I = el(e).offset : el(e).offset + el(e).ndof-1;

    M(I,I) = el(e).mass;

    for i = 1:el(e).ndof % .... Over local DoFs      
      for k = 1:el(e).nip % ... Over integration points

        fc = model.EvalConvFlux( el(e).sol(k) );
        R( I(i) ) = R( I(i) ) + el(e).qw(k) * el(e).dphi( i,k ) * fc;

      end
    end

  end


  % Ovel all faces
  for f = 1:Nf

    fc = NumConvFlux( fa(f).sol1, fa(f).sol2, model );

    % Left element
    e1    = fa(f).elnr1;
    I1    = el( e1 ).offset : el( e1 ).offset + el( e1 ).ndof-1;
    R(I1) = R(I1) - fa(f).phi1 * fc;

    % Right element
    e2    = fa(f).elnr2;
    I2    = el( e2 ).offset : el( e2 ).offset + el( e2 ).ndof-1;
    R(I2) = R(I2) + fa(f).phi2 * fc;

  end

  % \tilde{R} = M^{-1}R
  R = M\R;

end
