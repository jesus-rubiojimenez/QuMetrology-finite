% Mean square error for the estimation of two linear functions
% 
% The estimation scheme is a quantum sensing network with two qubits.
%
% Note that we use the trapezoidal rule 'trapz' for the inner parameter
% integrals because these have peaked integrands, while Simpson's Rule
% 'simps' is a better choice when this problem does not arise, which is
% the case for the outer parameter integrals.
%
% The 'simps' function employed in this code can be found in
%
%   Damien Garcia (2021). Simpson's rule for numerical integration 
%   (https://www.mathworks.com/matlabcentral/fileexchange/25754-simpson-s-rule-for-numerical-integration),
%   MATLAB Central File Exchange. Retrieved July 18, 2021.
clear

% Initial parameters
prior_mean1=pi/4;
prior_mean2=pi/4;
prior_width1=pi/2;
prior_width2=pi/2;
mu_max=1;

% Weighting matrix
WD=[1 0; 0 1]/2;

% Transformation representing the original parameters
%K=[1 0; 0 1];

% Transformation representing two linear functions
V=[2/sqrt(4+pi^2) 2/sqrt(5); pi/sqrt(4+pi^2) 1/sqrt(5)];
  
% Combination of linear transformation and weighting matrix
G=V*WD*V';

% Initial state
gamma_par=1; % Local strategy
%gamma_par=0; % Maximally entangled strategy
%gamma_par=0.530696; % Asymptotically optimal strategy
%gamma_par=0.3343605926149827; % Balanced startegy
initial_state=sparse([1 gamma_par gamma_par 1])'/sqrt(2+2*gamma_par^2);

% Generators
sigmaz=sparse([1 0; 0 -1]);
g1=kron(sigmaz,identity(2))/2;
g2=kron(identity(2),sigmaz)/2;

% Asymptotically optimal local POM (F = F_q, chapter 6)
proj1=sparse([-1 -1 1 1])'/2;
proj2=sparse([1 1 1 1])'/2;
proj3=sparse([1 -1 -1 1])'/2;
proj4=sparse([-1 1 -1 1])'/2;
proj_columns=[proj1';proj2';proj3';proj4']';

% Optimal single-shot POM (chapter 7)
% proj1=sparse([1i 1 1 -1i])'/2;
% proj2=sparse([-1i 1 1 1i])'/2;
% proj3=sparse([1i -1 1 1i])'/2;
% proj4=sparse([-1i -1 1 -1i])'/2;
% proj_columns=[proj1';proj2';proj3';proj4']';

% Parameter domain
dim_theta=100;
dim_theta_out=20;
a1=prior_mean1-prior_width1/2;
b1=prior_mean1+prior_width1/2;
theta1=linspace(a1,b1,dim_theta); % Inner parameter integrals
theta1_out=linspace(a1,b1,dim_theta_out); % Outer parameter integrals
a2=prior_mean2-prior_width2/2;
b2=prior_mean2+prior_width2/2;
theta2=linspace(a2,b2,dim_theta);
theta2_out=linspace(a2,b2,dim_theta_out);

% Monte Carlo sample size
tau_mc=200;

% State after encoding the parameters, final state and amplitudes
amplitudes=zeros(size(proj_columns,2),dim_theta,dim_theta);
amplitudes_sparse=zeros(dim_theta,dim_theta,size(proj_columns,2));
for z1=1:dim_theta
  for z2=1:dim_theta
    after_encoding=sparse(expm(-1i*(g1*theta1(z1)+g2*theta1(z2))))*initial_state;
    for x=1:size(proj_columns,2)
      povm_element=proj_columns(:,x);
      amplitudes_temp=sparse(povm_element)'*sparse(after_encoding);
      amplitudes(x,z1,z2)=amplitudes_temp;
      amplitudes_sparse(z1,z2,x)=amplitudes_temp;
    end
    
    % The second method of generating the amplitudes is included in
    % order to use sparse later in the code.
  end
end

% Likelihood function
likelihood=amplitudes.*conj(amplitudes);
if (1-sum(likelihood(:,1,1)))>1e-7
  error('The quantum probabilities do not sum to one.')
end
likelihood_sparse=amplitudes_sparse.*conj(amplitudes_sparse);
if (1-sum(likelihood_sparse(1,1,:),3))>1e-7
  error('The quantum probabilities (sparse version) do not sum to one.')
end

% Prior probability
prior=ones(dim_theta,dim_theta);
prior=prior/trapz(theta2,trapz(theta1,prior));
prior_out=ones(dim_theta_out,dim_theta_out);
prior_out=prior_out/trapz(theta2_out,trapz(theta1_out,prior_out));

% Bayesian mean square error
epsilon_out=zeros(dim_theta_out,dim_theta_out);
for index_out1=1:dim_theta_out
  for index_out2=1:dim_theta_out
    
    % Matching outer and inner parameter indices
    for y=1:dim_theta
      if theta1(y)>theta1_out(index_out1) || theta1(y)==theta1_out(index_out1)
        index_real1=y;
        break
      end
    end
    
    for z=1:dim_theta
      if theta2(z)>theta2_out(index_out2) || theta2(z)==theta2_out(index_out2)
        index_real2=z;
        break
      end
    end
    
    epsilon_n1=zeros(1,mu_max);
    epsilon_n2=zeros(1,mu_max);
    epsilon_n_offdia=zeros(1,mu_max);
    epsilon_n_sum=zeros(1,mu_max);
    for times=1:tau_mc

      % Prior density function
      prob_temp=sparse(prior);
      for runs=1:mu_max

        % (Monte Carlo) Outcome simulation
        prob_sim=likelihood(:,index_real1,index_real2);
        cumulative = cumsum(prob_sim); % Cumulative function
        prob_rand=rand; % Random selection
        auxiliar=cumulative-prob_rand;

        for x=1:size(proj_columns,2)
          if auxiliar(x)>0
            index_mc=x;
            break 
          end
        end

        % Posterior density function
        likesimulated=likelihood_sparse(:,:,index_mc);
        prob_temp=sparse(prob_temp.*likesimulated);
        normalisation=sparse(trapz(theta2,trapz(theta1,prob_temp,1),2));
        if normalisation>1e-16
          prob_temp=prob_temp/normalisation;
        else
          prob_temp=0;
        end
        prob_temp=sparse(prob_temp);

        % Bayes estimator for the first parameter
        theta_expe1=trapz(theta1,trapz(theta2,prob_temp,2).*theta1',1);
        theta2_expe1=trapz(theta1,trapz(theta2,prob_temp,2).*theta1'.^2,1);
        epsilon_n1(runs)=theta2_expe1-theta_expe1^2;

        % Bayes estimator for the second parameter
        theta_expe2=trapz(theta2,trapz(theta1,prob_temp,1).*theta2,2);
        theta2_expe2=trapz(theta2,trapz(theta1,prob_temp,1).*theta2.^2,2);
        epsilon_n2(runs)=theta2_expe2-theta_expe2^2;

        % Off-diagonal terms (the covariance matrix is symmetric)
        theta2_offdia=trapz(theta1,trapz(theta2,prob_temp.*theta2,2).*theta1',1);
        epsilon_n_offdia(runs)=theta2_offdia-theta_expe1*theta_expe2;

      end

      % Monte Carlo sum with transformation and weighting matrices
      epsilon_n_sum=epsilon_n_sum+G(1,1)*epsilon_n1+G(2,2)*epsilon_n2+2*G(1,2)*epsilon_n_offdia;
    end
    
    % Monte Carlo approximation
    epsilon_average=epsilon_n_sum/(tau_mc);
    for runs_out=1:mu_max
      epsilon_out(index_out1,index_out2,runs_out)=epsilon_average(runs_out);
    end
  end
end

% Outer integral
epsilon_trials=zeros(1,mu_max);
for runs_out=1:mu_max
  epsilon_temp=epsilon_out(:,:,runs_out);
  epsilon_trials(runs_out)=simps(theta2_out,simps(theta1_out,prior_out.*epsilon_temp));
end

% Observations
observations=1:1:mu_max;

% Fisher information matrix
F11=4*(initial_state'*g1^2*initial_state-(initial_state'*g1*initial_state)^2);
F12=4*(initial_state'*g1*g2*initial_state-(initial_state'*g1*initial_state)*(initial_state'*g2*initial_state));
F21=4*(initial_state'*g2*g1*initial_state-(initial_state'*g2*initial_state)*(initial_state'*g1*initial_state));
F22=4*(initial_state'*g2^2*initial_state-(initial_state'*g2*initial_state)^2);
F=[F11 F12; F21 F22];

% Quantum Cramer-Rao bound
qcrb=trace(G/F)./(observations);

% Save results
%save('qnetwork_results.txt','observations','epsilon_trials','qcrb','ascii')
