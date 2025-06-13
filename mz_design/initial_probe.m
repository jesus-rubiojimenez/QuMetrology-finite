function [initial_state] = initial_probe(state_sel)
% Common states in optical interferometry, where 'state_sel' is a number
% from 1 to 5 labelling the quantum probes
%
% (1) Coherent state: |alpha/sqrt(2),-i*alpha/sqrt(2)>
% (2) NOON state: (|N,0>+|0,N>)/sqrt(2)
% (3) Twin squeezed vacuum state: S_a(z)S_b(z)|0,0>
% (4) Squeezed entangled state: N (|z,0> + |0,z>)
% (5) Twin squeezed cat state: [N S(z)(|alpha>+|-alpha>)]\otimes2
%
% whose componentes are generated in the number basis of a Mach-Zehnder
% interferometer.
%
% The code is configured with the parameters
%
% (1) alpha=sqrt(2)
% (2) N=2
% (3) z=asinh(1)
% (4) z=log(2+sqrt(3))
% (5) alpha=0.960149, z=1.2145
%
% so that the nbar number of quanta that enters the interferometer is 2.
%
% The cutoff for the vectors are: (1) 20, (2) 2, (3) 50, (4) 60
% and (5) 50. These values are selected such that the numerical states
% are a reasonable approximation to the analytical kets.
  
% State parameters
nbar=2;
number=nbar; % Mean number of photons
alpha=sqrt(nbar); % Displacement parameter
zeta=asinh(sqrt(nbar/2)); % Squeezing parameter
zent=log(2+sqrt(3));

alphacat=0.960149; % Maximum Fisher information for the twin squeezed
zcat=1.2145;       % cat state

%alphacat=1.09048; % Same Fisher information for the twin squeezed cat
%zcat=1.1025;      % state and the squeezed entangled state

if state_sel==1
  num_cutoff=20; % Cutoff for states
elseif state_sel==2
  num_cutoff=number;
elseif state_sel==3 || state_sel==5
  num_cutoff=50;
elseif state_sel==4
  num_cutoff=60;
end
op_cutoff=num_cutoff+1; % Cutoff for operators

% Initial state
if state_sel==1
  initial_temp=sparse(displacement(op_cutoff,alpha)*vacuum(op_cutoff));
  initial_state=sparse(kron(initial_temp,vacuum(op_cutoff)));
  initial_state=beam_splitter(op_cutoff)*initial_state;
elseif state_sel==2
  initial_temp=sparse((creation(op_cutoff)^number)*vacuum(op_cutoff));
  initial_state=sparse(kron(initial_temp,vacuum(op_cutoff))+kron(vacuum(op_cutoff),initial_temp));
elseif state_sel==3
  initial_temp1=sparse(squeeze(op_cutoff,zeta)*vacuum(op_cutoff));
  initial_temp2=sparse(squeeze(op_cutoff,zeta)*vacuum(op_cutoff));
  initial_state=sparse(kron(initial_temp1,initial_temp2));
elseif state_sel==4
  initial_state=(kron(squeeze(op_cutoff,zent),identity(op_cutoff))+kron(identity(op_cutoff),squeeze(op_cutoff,zent)))*kron(vacuum(op_cutoff),vacuum(op_cutoff));
elseif state_sel==5
  initial_state=squeeze(op_cutoff,zcat)*(displacement(op_cutoff,alphacat)+displacement(op_cutoff,-alphacat))*vacuum(op_cutoff);
  initial_state=kron(initial_state,initial_state);
  initial_state=initial_state/sqrt((initial_state'*initial_state));
end
initial_state=initial_state/sqrt((initial_state'*initial_state)); end
