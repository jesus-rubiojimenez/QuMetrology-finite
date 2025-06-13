% Prior information analysis for single-parameter schemes
%
% This programme uses Bayes theorem to generate the posterior probability
%
% p(theta|m_1, ..., m_mu)
%
% for a flat prior and the likelihood function given by the Born rule.
% The initial state and the measurement scheme are those of a Mach-Zehnder
% interferometer, and they can be selected from the respective MATLAB
% functions in our interferometric toolbox by giving a value from 1 to 5
% for 'state_choice' and 'pom_choice'.
%
% Important observations:
%
% - The prior is defined over all the parameter domain, so that the
% symmetries of the likelihood that enable us to find the intrinsic
% width can be visualised.
%
% - The variables 'prior_mean_1shot' and 'prior_width_1shot' are needed
% to specify the optimal single-shot POM, but they do not affect the
% other measurement schemes.
%
% - The results for the prior information analysis in chapter 4 are
% recovered when we remove the extra phase shifts in the second option
% of our MATLAB function mz_pom(.) and we select it. This is because
% the results in chapter 4 were obtained for a prior between 0 and W,
% while the POMs included in our sample of codes are those associated
% with a prior centred around zero (see chapter 5).
clear

% State and POM options (see the respective codes in previous sections)
state_choice=2;
pom_choice=2;

% Initial state
initial_state=initial_probe(state_choice);

% Space cutoff
op_cutoff=sqrt(length(initial_state));

% Parameter domain
prior_mean=pi;
prior_width=2*pi; % Complete parameter domain
a=prior_mean-prior_width/2;
b=prior_mean+prior_width/2;
dim_theta=1000;
theta=linspace(a,b,dim_theta);

% Simulation of the unknown true value
index_real=160;

% Measurement scheme
prior_mean_1shot=0;
prior_width_1shot=pi/2;
[outcomes_space,proj_columns] = mz_pom(state_choice,pom_choice,prior_width_1shot,prior_mean_1shot);

% State after the phase shift, final state and amplitudes
amplitudes=zeros(length(outcomes_space),dim_theta);
for z=1:dim_theta
  after_phase_shift=sparse(phase_shift_diff(op_cutoff,theta(z))*initial_state);
  for x=1:length(outcomes_space)
    pom_element=proj_columns(:,x);
    amplitudes(x,z)=sparse(pom_element'*after_phase_shift);
  end
end

% Likelihood function (using the Born rule)
likelihood=amplitudes.*conj(amplitudes);

% Prior density function
prior=ones(1,dim_theta);
prior=prior/trapz(theta,prior);

% Updating via Bayes theorem
prob_temp=prior;
for runs=1:100
  
  % Simulation of an interferometric experiment
  prob_sim=likelihood(:,index_real);
  cumulative = cumsum(prob_sim); % Cumulative function
  prob_rand=rand; % Random selection
  auxiliar=cumulative-prob_rand;

  for x=1:length(outcomes_space)
    
    if auxiliar(x)>0
      index=x;
      break
    end
  end

  % Posterior density function
  prob_temp=sparse(prob_temp.*likelihood(index,:));
  if trapz(theta,prob_temp)>1e-16
    prob_temp=prob_temp./trapz(theta,prob_temp);
  else
  prob_temp=0;
  end

  % Posterior probability plots
  if runs==1
    plot(theta,prob_temp,'k-','LineWidth',2.5)
    hold on
  elseif runs==2; plot(theta,prob_temp,'k-','LineWidth',2.5)
  elseif runs==10; plot(theta,prob_temp,'k-','LineWidth',2.5)
  elseif runs==100; plot(theta,prob_temp,'k-','LineWidth',2.5)
    hold off
  end
end

% Plot specifications
grid
fontsize=21;
set(gcf,'units','points','position',[250,50,550,400])
xlabel('$\theta$','Interpreter','latex','FontSize',fontsize)
ylabel('$p(\theta | \textbf{\textit{m}})$','Interpreter','latex','FontSize',fontsize)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlim([min(theta) max(theta)])
set(gca,'FontSize', fontsize,'FontName','Times New Roman')
