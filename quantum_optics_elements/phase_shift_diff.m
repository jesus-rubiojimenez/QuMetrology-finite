function [utheta] = phase_shift_diff(dimension, theta)
% Matrix representation of the unitary encoding of the unknown parameter
% ’theta’ (difference of phase shifts).
utheta=expm(-1i*theta*j3schwinger(dimension));
utheta=sparse(utheta); end
