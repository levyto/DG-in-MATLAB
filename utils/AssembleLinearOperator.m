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
% Description:  Assemble the linear operator L in the semi-discretized problem %
%                                                                              %
%                 dU                                                           %
%                ---- = R(U) = L*U                                             %
%                 dt                                                           %
%                                                                              %
%               i.e. we are implicitly assuming here that we are solving       %
%               linear PDE, i.e. linear advection equation.                    %
%                                                                              %
%                                                                              %
% Input:       el ... Array of elements                                        %
%              fa ... Array of faces with                                      %
%           model ... Mathematical model - Here ScalarLinearAdvection1D        %
%                                                                              %
% Output:       L ... Matrix of the linear operator                            %
%               M ... Global mass matrix                                       %
%               D ... Global differentiation matrix                            %
%               F ... Global flux matrix                                       %
%                                                                              %
% ---------------------------------------------------------------------------- %
function [L, M, D, F] = AssembleLinearOperator(el, fa, model)

  Ne = el(1).Ne;
  Nf = fa(1).Nf;
  DoF = el(1).DoF;

  % Global operator matrix
  L = zeros(DoF);

  % Global mass matrix
  M = zeros(DoF);

  % Global differentiation matrix
  D = zeros(DoF);

  % Global flux matrix
  F = zeros(DoF);
  % FL = zeros(DoF);
  % FR = zeros(DoF);

  % Ovel all elements
  for e = 1:Ne

    % Indices corresponding to element e in the global matrix
    I = el(e).offset : el(e).offset + el(e).ndof-1;

    % Local mass matrix
    M(I,I) = el(e).mass;

    % Local differentiation matrix
    for i = 1:el(e).ndof  %...... Over local DoFs  
      for j = 1:el(e).ndof % .... Over local DoFs      
        for k = 1:el(e).nip % ... Over integration points
          D(I(i),I(j)) = D(I(i),I(j)) + el(e).qw(k) * el(e).dphi( i,k ) * el(e).phi( j,k ); 
        end
      end
    end

  end

  D = D*model.a;


  % Ovel all faces
  for f = 1:Nf

    % Indices corresponding to the upwind flux
    E = fa(f).elnr1;
    J = el(e).offset : el(e).offset + el(e).ndof-1;

    % Flux on the left (i.e. right interface of elnr2, hence the plus)
    % ----------------
    e = fa(f).elnr1;
    % Indices corresponding to element e in the global matrix 
    I = el(e).offset : el(e).offset + el(e).ndof-1; % ... (same as J in fact)
    % Left part of the ocal flux matrix
    for i = 1:el(e).ndof  %...... Over local DoFs  
      for j = 1:el(e).ndof % .... Over local DoFs 
        F(I(i),J(j)) = F(I(i),J(j)) + fa(f).phi1(i) * fa(f).phi1(j);
        % FL(I(i),J(j)) = fa(f).phi1(i) * fa(f).phi1(j);
      end
    end

    % Flux on the right (i.e. left interface of elnr2, hence the minus)
    % -----------------
    e = fa(f).elnr2;
    % Indices corresponding to element e in the global matrix
    I = el(e).offset : el(e).offset + el(e).ndof-1;
    % Left part of the ocal flux matrix
    for i = 1:el(e).ndof  %...... Over local DoFs  
      for j = 1:el(e).ndof % .... Over local DoFs 
        F(I(i),J(j)) = F(I(i),J(j)) - fa(f).phi2(i) * fa(f).phi1(j);
        % FR(I(i),J(j)) = - fa(f).phi2(i) * fa(f).phi1(j);
      end
    end

  end

  F = F*model.a;
  

  L = - M \ ( -D + F );

end