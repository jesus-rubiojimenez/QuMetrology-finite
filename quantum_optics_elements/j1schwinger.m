function [j1] = j1schwinger(dimension)
% Matrix representation of the J1 operator (Jordan-Schwinger map).
j1=0.5*(kron(creation(dimension),creation(dimension)') + kron(creation(dimension)',creation(dimension)));
j1=sparse(j1); end
