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
% Description:  ?????????????????????????????????????????????????????????????  %
%               ?????????????????????????????????????????????????????????????  %
%               ?????????????????????????????????????????????????????????????  %
%                                                                              %
% ---------------------------------------------------------------------------- %
function ut_function()

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%                               BEGIN                                %%%%%
  correct = true;
  name = "Title of the unittest1"

  % PREPARATION:
  #
  #
  #


  % TEST:
  if #COND#
    correct = false;
  end

  msg(correct, name, dbstack);
  %%%%%                                END                                 %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%                               BEGIN                                %%%%%
  correct = true;
  name = "Title of the unittest2"

  % PREPARATION:
  #
  #
  #


  % TEST
  if #COND#
    correct = false;
  end

  msg(correct, name);
  %%%%%                                END                                 %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end