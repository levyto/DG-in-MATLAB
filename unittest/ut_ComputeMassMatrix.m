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
% Description:  Unit test for computing mass matrix on one element. When       %
%               using orthogonal basis functions, the exact matrix is          %
%               diagonal.                                                      %
%                                                                              %
% ---------------------------------------------------------------------------- %
function ut_ComputeMassMatrix()

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%                               BEGIN                                %%%%%
  correct = true;
  name = "Mass matrix with LG integration is diagonal";

  % PREPARATION ***************************************************************
  
  order = 50;
  order_int = 2*order;

  [x,w] = GaussLegendreRule( order_int, -1, 1 );
  P = LegendrePolynomial( x, order, -1, 1 );
  M = ComputeMassMatrix(P,w);

  MmDiag = M - diag(M);
  err = sum(sum(MmDiag));


  % TEST:
  if err > 1e-10;
    correct = false;
  end

  msg(correct, name, dbstack);
  %%%%%                                END                                 %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%                               BEGIN                                %%%%%
  correct = true;
  name = "Row-sum with LG integration is exact";

  % PREPARATION:

  % Diagonal elements of the exact mass matrix for p = 50 on [-1,1]
  dof = order + 1;
  m = zeros(dof,1);
  for i = 1:dof
    l = i-1;
    m(i) = 2/(2*l + 1);
  end

  % Compute diag of lumped mass matrix
  mM  = sum(abs(M),2);
  err = sum(abs(m - mM));


  %%% TEST: 
  if err > 1e-10;
    correct = false;
  end

  msg(correct, name);
  %%%%%                                END                                 %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%                               BEGIN                                %%%%%
  correct = true;
  name = "Mass matrix with LGL integration is diagonal";

  % PREPARATION ***************************************************************

  [x,w] = LobattoGaussLegendreRule( order_int, -1, 1 );
  P = LegendrePolynomial( x, order, -1, 1 );
  M = ComputeMassMatrix(P,w);

  MmDiag = M - diag(M);
  err = sum(sum(MmDiag));


  % TEST:
  if err > 1e-10;
    correct = false;
  end

  msg(correct, name);
  %%%%%                                END                                 %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%                               BEGIN                                %%%%%
  correct = true;
  name = "Row-sum with LGL integration is exact";

  % PREPARATION:

  % Compute diag of lumped mass matrix
  mM  = sum(abs(M),2);
  err = sum(abs(m - mM));


  %%% TEST: 
  if err > 1e-10;
    correct = false;
  end

  msg(correct, name);
  %%%%%                                END                                 %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end