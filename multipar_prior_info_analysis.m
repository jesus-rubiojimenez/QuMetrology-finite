% Two-parameter prior information analysis
%
% Prior information analysis for a qubit sensing network. The basic logic
% of the method parallels that for the single-parameter case (see appendix
% B.5 and chapter 6 for more details).
  
% Initial parameters
prior_mean1=pi;
prior_mean2=prior_mean1;
prior_width1=2*pi;
prior_width2=2*pi;
mu_max=100;

% True values for the unknown parameters
theta1_real=1;
theta2_real=2;

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
dim_theta=200;
a1=prior_mean1-prior_width1/2;
b1=prior_mean1+prior_width1/2;
theta1=linspace(a1,b1,dim_theta);
a2=prior_mean2-prior_width2/2;
b2=prior_mean2+prior_width2/2;
theta2=linspace(a2,b2,dim_theta);

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
  error('The quantum probabilities do not sum to one.')
end

% Prior probability
prior=ones(dim_theta,dim_theta);
prior=prior/trapz(theta2,trapz(theta1,prior));

% Simulation of the true values for the unkonwn parameters
for y=1:dim_theta
  if theta1(y)>theta1_real || theta1(y)==theta1_real
    index_real1=y;
    break
  end
end

for y=1:dim_theta
  if theta2(y)>theta2_real || theta2(y)==theta2_real
    index_real2=y;
    break
  end
end

% Bayesian simulation
outcomes=zeros(1,mu_max);
for runs=1:mu_max

  % Simulation of the experimental outputs
  prob_sim=likelihood(:,index_real1,index_real2);
  cumulative1 = cumsum(prob_sim); % Cumulative function
  prob_rand=rand; % Random selection
  auxiliar=cumulative1-prob_rand;

  for x=1:size(proj_columns,2)
    if auxiliar(x)>0
      index=x;
      break
    end
  end
    
  outcomes(runs)=index;
end

% Prior density function
prob_temp=prior;
for runs=1:mu_max

  % Updated posterior density function
  ytemp=outcomes(runs);

  likesimulated=likelihood_sparse(:,:,ytemp);
  prob_temp=sparse(prob_temp.*likesimulated);
  prob_norm=sparse(trapz(theta2,trapz(theta1,prob_temp,1),2));
  if prob_norm>1e-16
    prob_temp=prob_temp/prob_norm;
  else
    prob_temp=0;
  end
  prob_temp=sparse(prob_temp);
end

% Plot of the posterior
contour(theta1',theta2',prob_temp,'LevelStep',0.1,'Fill','on')
xticks([0 pi/4 pi/2 3*pi/4 pi 5*pi/4 3*pi/2 7*pi/4 2*pi])
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4', '3\pi/2','7\pi/4','2\pi'})
yticks([0 pi/4 pi/2 3*pi/4 pi 5*pi/4 3*pi/2 7*pi/4 2*pi])
yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4', '3\pi/2','7\pi/4','2\pi'})
xt = get(gca, 'XTick');
fontsize=32;
set(gca,'FontSize', fontsize,'FontName','Times New Roman');
yt = get(gca, 'YTick');
set(gca,'FontSize', fontsize,'FontName','Times New Roman');
grid
