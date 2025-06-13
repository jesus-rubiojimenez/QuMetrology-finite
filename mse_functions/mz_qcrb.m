function [qcrb] = mz_qcrb(initial_state,mu_max)
% Quantum Cramer-Rao bound as a function of the number of trials, where
% 'initial_state' is a pure state for the Mach-Zehnder interferometer
% and 'mu_max' is the maximum number of repetitions.
  
% Number of repetitions
observations=1:1:mu_max;

% Space cutoff (for a single mode)
op_cutoff=sqrt(length(initial_state));

% Quantum Fisher information (pure states and unitary encoding)
expectation_n=initial_state'*j3schwinger(op_cutoff)*initial_state;
expectation_n2=initial_state'*j3schwinger(op_cutoff)^2*initial_state;
qfi=4*(expectation_n2-expectation_n^2);

% Do we have information?
if qfi==0
  disp('The Quantum Fisher information is zero.')
  return
end

% Quantum Cramer-Rao Bound
qcrb=1./(observations*qfi);
end
