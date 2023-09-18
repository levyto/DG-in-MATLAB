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
% Description:  Evaluate numerical flux based on left and right state          %
%                                                                              %
% Input:       w1 ... Solution evaluated on the left side of the interface     %
%              w2 ... Solution evaluated on the right side of the interface    %
%           model ... Mathematical model defining the flux                     %
%                                                                              %
% Output:    flux ... Value of the numerical flux at the interface             %
%                                                                              %
% ---------------------------------------------------------------------------- %
function flux = NumConvFlux(w1, w2, model)

  f1 = model.EvalConvFlux(w1);
  f2 = model.EvalConvFlux(w2);

  % Lax-Friedrichs
  lambda = max( model.EvalMaxEigenValue(w1), model.EvalMaxEigenValue(w2) );
  flux = 0.5*( f1 + f2 ) + 0.5*lambda*( w1 - w2 );

end