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
% Description:  Run all unit tests                                             %
%                                                                              %
% ---------------------------------------------------------------------------- %
function RunAllTestScripts()

  addpath("../src");

  fprintf("\n%s\n",repmat("=",1,80));
  fprintf(" Running unit tests:"); 
  tic
  fprintf("\n%s\n",repmat("=",1,80));

  % ****************************************************************************


  ut_LegendrePolynomial();
  ut_GaussLegendreRule();
  ut_LobattoGaussLegendreRule();
  ut_ComputeMassMatrix();


  % ****************************************************************************

  fprintf("\n%s\n",repmat("=",1,80));
  fprintf(" Running unit tests finished.\n Execution time: "); 
  toc
  fprintf("%s\n\n",repmat("=",1,80));

  return

end