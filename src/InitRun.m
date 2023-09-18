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
% Description:  Initialize the solution coefficient vector,                    %
%               perform L2-projection of the initial condition,                %
%               and reconstruct the solution on each element and face.         %
%               Display info about the size of the problem to be solved.       %
%                                                                              %
% Input:       el ... Array of elements                                        %
%              fa ... Array of faces                                           %
%             sol ... Global coefficient vector to be computed                 %
%           input ... Input parameters object                                  %
%                                                                              %
% Output:     sol ... Global coefficient vector based on initial condition     %
%              el ... Array of elements with stored initial condition          %
%              fa ... Array of faces with stored initial condition             %
%                                                                              %
% ---------------------------------------------------------------------------- %
function [sol, el, fa] = InitRun(el, fa, input)

  % Initialize size of coefficient vector
  sol = zeros(el(1).DoF,1);

  % L2-Projection of initial condition
  sol = InitialConditionL2Projection(el, fa, sol, input.ic);

  % Reconstruct the solution on each element based on the coefficients
  [el, fa] = ReconstructSolution(el, fa, sol);
  
  % Print info
  fprintf("\n")

  fprintf("%-25s%d\n","  Total DoFs:",el(1).DoF);
  fprintf("%-25s%d\n","  Total int. points:",el(1).NIP);
  fprintf("%-25s%d\n","  Number of faces:",fa(1).Nf);

  fprintf("\n")

  % Plot initial condition
  if input.showic == 1

    PlotSolution(el, fa, sol);
    fprintf("  Paused after plotting the initial condition.\n");
    fprintf("  Waiting for input ...\n");

    pause

  end

end