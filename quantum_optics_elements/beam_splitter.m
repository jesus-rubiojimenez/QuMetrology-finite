function [v] = beam_splitter(dimension)
% Matrix representation of a 50:50 beam splitter.
v=expm(-1i*0.5*pi*sparse(j1schwinger(dimension)));
v=sparse(v); end
