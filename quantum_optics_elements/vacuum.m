function [zero] = vacuum(dimension)
% Vacuum state for a single mode.
temp=identity(dimension);
zero=temp(:,1);
zero=sparse(zero); end
