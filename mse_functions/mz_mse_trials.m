function[epsilon_trials]=mz_mse_trials(state_choice,pom_choice,prior_width,prior_mean,mu_max)
% Bayesian mean square error
%
% This programme calculates the mean square error as a function of the
% number of repetitions.
%
% To run it, we need to specify the variables 'state_choice', which
% labels the initial state of a Mach-Zehnder interferometer; 'pom_choice',
% which selects the measurement scheme; 'phase_width', which is the
% width of a flat prior probability; 'phase_mean', which is the centre
% of its domain; and 'mu_max', which is the maximum number of trials.
%
% Note that this code relies on other MATLAB functions of our numerical
% toolbox. The algorithm in this section has been exploited to calculate
% the mean square error for all the single-parameter cases treated in this
% thesis, including the ideal schemes for optical interferometry studied
% in chapters 4 and 5, the calculation of the Taylor error to verify the
% validity of our squared approximation in appendix A, our lossy analysis
% in chapter 8 and our analysis of the elapsed time, also in chapter 8.
  
% Seed for the random generator
rng('shuffle')
  
% Initial state
initial_state=initial_probe(state_choice);

% Space cutoff
op_cutoff=sqrt(length(initial_state));

% Parameter domain
a=prior_mean-prior_width/2;
b=prior_mean+prior_width/2;
dim_theta=1250;
theta=linspace(a,b,dim_theta);
num_steps=125;
step=round(dim_theta/num_steps);
if step-round(step)~=0
  disp('Error: dim_theta divided by num_steps must be an integer.')
  return
elseif num_steps<3
  disp('Error: the approximation for the external theta integral needs three rectangles at least.')
return
end

% Monte Carlo sample size
tau_mc=1250;

% Measurement scheme
[outcomes_space,proj_columns] = mz_pom(state_choice,pom_choice,prior_width, prior_mean);

% State after the phase shift, final state and amplitudes
amplitudes=zeros(length(outcomes_space),dim_theta);
for z=1:dim_theta
  after_phase_shift=sparse(phase_shift_diff(op_cutoff,theta(z))*initial_state);
  for x=1:length(outcomes_space)
    pom_element=proj_columns(:,x);
    amplitudes(x,z)=sparse(pom_element'*after_phase_shift);
  end
end

% Likelihood function
likelihood=amplitudes.*conj(amplitudes);
disp('The likelihood function has been created.')
if (1-sum(likelihood(:,1)))>1e-7
  error('The quantum probabilities do not sum to one.')
end

% Prior probability
prior=ones(1,dim_theta);
prior=prior/trapz(theta,prior);

% Bayesian inference
epsilon_bar=0;
for index_real=1:step:dim_theta
  epsilon_n=zeros(1,mu_max); % Preallocate vector
  epsilon_n_sum=zeros(1,mu_max);
  for times=1:tau_mc

    % Prior density function
    prob_temp=prior;
    for runs=1:mu_max
  
      % (Monte Carlo) Interferometric simulation
      prob_sim1=likelihood(:,index_real);
      cumulative1 = cumsum(prob_sim1); % Cumulative function
      prob_rand1=rand; % Random selection
      auxiliar1=cumulative1-prob_rand1;
    
      for x=1:length(outcomes_space)
        if auxiliar1(x)>0
          index1=x;
          break
        end
      end

      % Posterior density function
      prob_temp=sparse(prob_temp.*likelihood(index1,:));
      if trapz(theta,prob_temp)>1e-16
        prob_temp=prob_temp./trapz(theta,prob_temp);
      else
        prob_temp=0;
      end

      % Experimental square error
      theta_expe=trapz(theta,prob_temp.*theta);
      theta2_expe=trapz(theta,prob_temp.*theta.^2);
      epsilon_n(runs)=theta2_expe-theta_expe^2;
    end

    % Monte Carlo sum
    epsilon_n_sum=epsilon_n_sum+epsilon_n;
  end

  % Monte Carlo approximation for the Bayesian error
  epsilon_average=epsilon_n_sum/(tau_mc);
  epsilon_bar=epsilon_bar+epsilon_average*prior(index_real)*(theta(2*step)-theta(step));
end
epsilon_trials=epsilon_bar;
