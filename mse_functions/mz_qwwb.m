function [qwwb] = mz_qwwb(initial_state,phase_width,mu_max)
% Quantum Weiss-Weinstein bound as a function of the number of
% trials, where 'initial_state' is a pure state for the Mach-Zehnder,
% interferometer 'phase_width' is the width of the parameter domain
% and ’mu_max’ is the maximum number of repetitions.
  
% Space cutoff (for a single mode)
op_cutoff=sqrt(length(initial_state));

% Parameter domain
W=phase_width;
dim_theta=1000;
h=linspace(0,W,dim_theta);

% Quantum Weiss-Weinstein Bound
find_supremum=zeros(mu_max,length(h));
qwwb=zeros(1,mu_max);
for runs=1:mu_max
  for z=1:length(h)
    zeta=initial_state'*phase_shift_diff(op_cutoff,h(z))*initial_state;
    zeta2=initial_state'*phase_shift_diff(op_cutoff,2*h(z))*initial_state;
    fid_function=abs(zeta)^2;
    find_supremum(runs,z)=h(z)^2*(1-h(z)/W)^2*fid_function^(2*runs)/(2*fid_function^runs-2*(1-2*h(z)/W)*real(zeta^(2*runs)*conj(zeta2)^runs));
  end
  qwwb(runs)=max(find_supremum(runs,:));
end
end
