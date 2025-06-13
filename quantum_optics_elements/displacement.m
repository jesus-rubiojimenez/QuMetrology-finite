function [displ] = displacement(dimension,alpha)
% Matrix representation of the displacement operator, where 'alpha' is
% the amount of displacement.
displ=expm(alpha*creation(dimension)-conj(alpha)*creation(dimension)');
displ=sparse(displ); end
