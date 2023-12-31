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
% Description:  Unittests for computing Legendre polynomials.                  %
%               Check correctness for the first 5 basis functions.             %
%                                                                              %
% ---------------------------------------------------------------------------- %
function ut_LegendrePolynomial()

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%                               BEGIN                                %%%%%
  correct = true;
  name = "Pointwise exactness of functions";

  % PREPARATION:

  x = -1:0.1:1;
  N = 4;

  P = zeros(N+1,length(x));
  P(1,:) = 1;
  P(2,:) = x;
  P(3,:) = 1/2*(3*x.^2 - 1);
  P(4,:) = 1/2*(5*x.^3 - 3*x);
  P(5,:) = 1/8*(35*x.^4 - 30*x.^2 + 3);

  P1 = LegendrePolynomial(x,N,-1,1);

  f = sum(sum(abs(P-P1)));


  % TEST:
  if ( f > 1e-13)
    correct = false;
  end

  msg(correct, name, dbstack);
  %%%%%                                END                                 %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%                               BEGIN                                %%%%%
  correct = true;
  name = "Pointwise exactness of derivatives";

  % PREPARATION:

  dP = zeros(N+1,length(x));
  dP(1,:) = 0;
  dP(2,:) = 1;
  dP(3,:) = 3*x;
  dP(4,:) = 1/2*(15*x.^2 - 3);
  dP(5,:) = 1/8*(140*x.^3 - 60*x);

  [~, dP1] = LegendrePolynomial(x,N,-1,1);

  df = sum(sum(abs(P-P1)));


  % TEST:
  if (df > 1e-13)
    correct = false;
  end

  msg(correct, name);
  %%%%%                                END                                 %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end