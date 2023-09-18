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
% Description:  Compute Legendre basis functions of order N and its            %
%               derivatives at point/points x defined on interval [-1,1].      %
%               Computation is done by the Bonnet's recursive formula.         %
%                                                                              %
% Input:        x ... Single point or vector of points where the functions     %
%                     have to be evaluated. x is in [a,b] and transformed      %
%                     to [-1,1] before the computation                         %
%               N ... Highest order of the polynomials                         %
%                                                                              %
% Output:       P ... Legendre pasis functions evaluated at point/points x,    %
%                     P and x have the same length                             %
%              dP ... Derivatives of Legendre pasis functions evaluated at     %
%                     point/points x, dP and x have the same length            %
%                                                                              %
% ---------------------------------------------------------------------------- %
function [P, dP] = LegendrePolynomial(x, N, a, b)

  % Funtions are defined on [-1, 1]: Transform points from [a,b] to [-1,1]
  if nargin > 2
    x = (2*x - (b+a))/(b-a);
  end

   P = zeros(N+1,length(x));
  dP = zeros(N+1,length(x));

  % n = 0
   P(1,:) = 1;
  dP(1,:) = 0;
  if N == 0; return; end

  % n = 1
   P(2,:) = x;
  dP(2,:) = 1;

  % n > 1
  n = 1;
  while n <= N-1

    % indices of Legendre polynomials
    i = [n+1, n, n-1];
    % indexing from 1
    i = i+1;

     P(i(1),:) = (2*n+1)/(n+1) * x .*  P(i(2),:) - n/(n+1) *  P(i(3),:);
    dP(i(1),:) = (2*n+1)/(n+1) * P(i(2),:) ...
               + (2*n+1)/(n+1) * x .* dP(i(2),:) - n/(n+1) * dP(i(3),:);

    n = n+1;

  end

  if nargin > 2
    % Scaling due to the transformation from [a,b] to [-1,1]
    dP(:,:) = dP * 2/(b-a);
  end

end
