function [Sopt,sopt,soptvec_columns,bayes_bound,rho,pk,psik,rhobar,rhobarnew] = mz_optimal_1trial(initial_state,phase_width,phase_mean)
% Optimal single-shot strategy, where 'initial_state' is a pure state
% for the Mach-Zehnder interferometer, 'phase_width' is the width of the
% parameter domain and 'phase_mean' is its centre.
%
% This programme calculates:
%
%   a) the optimal quantum estimator 'Sopt'
%   b) the estimates 'sopt' for the unknown parameter given by the
%      spectrum of 'Sopt'
%   c) the optimal projective measurement for a single trial given by
%      the eigenvectors 'soptvec_columns' of 'Sopt'
%   d) the optimal single-shot mean square error 'bayes_bound'
%   e) the zero-th quantum moment of the transformed density matrix 'rho'
%   f) 'rho' in its diagonal basis, denoted by 'pk'
%   g) the matrix 'psik' whose columns are the eigenvectors of 'rho'
%   h) the first quantum moment of the transformed density matrix 'rhobar'
%   i) 'rhobar' in the eigenbasis of 'rho', denoted by 'rhobarnew'

% Calculation of 'rho' and 'rhobar'
index=1;
kvec=zeros(1,length(initial_state)^2);
lvec=zeros(1,length(initial_state)^2);
for x1=1:sqrt(length(initial_state))
  for y1=1:sqrt(length(initial_state))
    for z1=1:sqrt(length(initial_state))
      for t1=1:sqrt(length(initial_state))
        if (x1-1)-(y1-1)+(t1-1)-(z1-1)==0
          K=phase_width;
          L=phase_mean*phase_width;
        else
          comp_temp=(x1-1)-(y1-1)+(t1-1)-(z1-1);
          exp_temp=exp(-1i*comp_temp*phase_mean/2);
          sin_temp=sin(comp_temp*phase_width/4);
          cos_temp=cos(comp_temp*phase_width/4);
          K=4*exp_temp*sin_temp/comp_temp;
          L=exp_temp*(4*phase_mean*sin_temp/comp_temp+1i*2*phase_width*cos_temp/comp_temp - 1i*8*sin_temp/comp_temp^2);
        end
        kvec(index)=K/phase_width;
        lvec(index)=L/phase_width;
        index=index+1;
      end
    end
   end
end

kmat=sparse(vec2mat(kvec,sqrt(length(kvec))));
lmat=sparse(vec2mat(lvec,sqrt(length(lvec))));
initial_rho=kron(initial_state,initial_state');
rho=initial_rho.*kmat; rhobar=initial_rho.*lmat;
rho=full(rho); rhobar=full(rhobar);

% Eigenvalues and eigenvectors of 'rho'
[psik, pk] = eigs(rho,rank(rho));
psik=sparse(psik);
pk=sparse(pk);

pkvec=zeros(1,length(pk));
for x=1:length(pk)
  pkvec(x)=pk(x,x);
end

% 'rhobar' in the eigenbasis of 'rho'
rhobarnew=psik'*rhobar*psik;

% Optimal single-shot strategy: projectors and outcomes
Sopt_temp=zeros(length(pkvec),length(pkvec));
for a=1:length(pkvec)
  for b=1:length(pkvec)
    if pkvec(a)+ pkvec(b)>0
      Sopt_temp(a,b)=2*rhobarnew(a,b)/(pkvec(a)+pkvec(b));
    end
  end
end

Sopt_temp=sparse(Sopt_temp);
Sopt=psik*Sopt_temp*psik';
Sopt=full(Sopt);

[soptvec_columns, sopt_temp] = eigs(Sopt,rank(Sopt));
sopt=zeros(1,length(sopt_temp));
for x=1:length(sopt_temp)
  sopt(x)=sopt_temp(x,x);
end
soptvec_columns=sparse(soptvec_columns);
sopt=sparse(sopt);
if imag(sopt)<1e-5; sopt=real(sopt);
else
  error('The estimates of the unknown parameter must be real. Check the cutoff in the intermediate calculations.')
return
end

% Phase domain
phase=linspace(phase_mean-phase_width/2,phase_mean+phase_width/2,1000);

% Optimal single-shot mean square error
bayes_bound=trapz(phase,phase.*phase)/phase_width - trace(Sopt*Sopt*rho);
if imag(bayes_bound)<1e-10; bayes_bound=real(bayes_bound);
else
  error('The mean square error must be real. Check the cutoff in the intermediate calculations.')
  return
end
end
