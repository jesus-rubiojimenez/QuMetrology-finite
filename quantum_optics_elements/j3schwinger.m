function [j3] = j3schwinger(dimension)
% Matrix representation of the J3 operator (Jordan-Schwinger map).
j3=0.5*(kron(creation(dimension)*creation(dimension)',identity(dimension))-kron(identity(dimension),creation(dimension)*creation(dimension)'));
j3=sparse(j3); end
