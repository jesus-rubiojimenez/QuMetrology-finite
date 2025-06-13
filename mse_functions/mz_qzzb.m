function [qzzb] = mz_qzzb(initial_state,phase_width,mu_max)
% Quantum Ziv-Zakai bound as a function of the number of trials, where
% 'initial_state' is a pure state for the Mach-Zehnder interferometer,
% 'phase_width' is the width of the parameter domain and 'mu_max' is
% the maximum number of repetitions.

% Space cutoff (for a single mode)
op_cutoff=sqrt(length(initial_state));

% Parameter domain
W=phase_width;
dim_theta=1000;
theta=linspace(0,W,dim_theta);

% Fidelity
fidelity=zeros(dim_theta,1);
for z=1:dim_theta
  after_phase_shift=sparse(phase_shift_diff(op_cutoff,theta(z))*initial_state);
  fidelity(z)=abs(initial_state'*after_phase_shift)^2;
end

% Quantum Ziv-Zakai Bound integrand
integrand=zeros(dim_theta,mu_max);
for runs=1:mu_max
  for z=1:dim_theta
    integrand(z,runs)=0.5.*theta(z).*(1-theta(z)./W).*(1-sqrt(1-fidelity(z).^runs));
  end
end
integrand=sparse(integrand);

% Quantum Ziv-Zakai Bound
qzzb=trapz(theta,integrand,1);
end
