function [outcomes,proj_columns] = mz_pom(state_choice,pom_choice,phase_width,phase_mean)
% Outcomes and POM elements of five projective measurement schemes:
%
%   1) Optimal single-shot POM
%   2) 50:50 beam splitter + photon counting
%   3) 50:50 beam splitter + measurement of quadratures rotated by pi/8
%   4) Undoing the preparation of the initial state + photon counting
%   5) 50:50 beam splitter + parity measurements
%
% where 'state_choice' labels the initial state, 'pom_choice' selects one
% of the previous measurement schemes, 'phase_width' is the width of the
% phase domain and 'phase_mean' is its centre.
%
% Some extra phase shifts that are assumed to be known have been added
% to 2) - 5) in order to make the strategy optimal when the prior is
% centred around zero.
  
% Space cutoff (for a single mode)
op_cutoff=sqrt(length(initial_probe(state_choice)));

if pom_choice==1
  % 1) Optimal single-shot POM
  [~,outcomes,proj_columns,~,~,~,~,~,~]=mz_optimal_1trial(initial_probe(state_choice),phase_width,phase_mean);

elseif pom_choice==2
  % 2) 50:50 beam splitter + photon counting
  
  % Observable quantity (number of photons at each port)
  observable=kron(creation(op_cutoff)*creation(op_cutoff)',creation(op_cutoff)*creation(op_cutoff)');
  [proj_columns,outcomes_temp]=eig(full(observable));
  outcomes=zeros(1,length(outcomes_temp));

  for x=1:length(outcomes_temp)
    outcomes(x)=outcomes_temp(x,x);
   end

  % Extra phase shift
  odd_shift=kron(identity(op_cutoff),expm(1i*(pi/2)*creation(op_cutoff)*creation(op_cutoff)'));
  even_shift=kron(identity(op_cutoff),expm(1i*(pi/4)*creation(op_cutoff)*creation(op_cutoff)'));

  if state_choice==1
    optimal_shift=odd_shift;
  else
  optimal_shift=even_shift;
  end
  
  % Effect of the 50:50 beam splitter
  proj_columns=optimal_shift'*beam_splitter(op_cutoff)*proj_columns;

elseif pom_choice==3
  % 3) 50:50 beam splitter + measurement of quadratures rotated by pi/8
  
  % Observable quantity
  if state_choice==1
    error('The quadrature POM is not available for coherent states.')
  else
    phasequad1=pi/8;
    phasequad2=phasequad1;
  end
  quad1=(creation(op_cutoff)*exp(1i*phasequad1)+creation(op_cutoff)'*exp(-1i*phasequad1))/sqrt(2);
  quad2=(creation(op_cutoff)*exp(1i*phasequad2)+creation(op_cutoff)'*exp(-1i*phasequad2))/sqrt(2);
  observable=kron(quad1,quad2);
  [proj_columns,outcomes_temp]=eig(full(observable));
  outcomes=zeros(1,length(outcomes_temp));
  for x=1:length(outcomes_temp)
    outcomes(x)=outcomes_temp(x,x);
  end
  
  % Extra phase shift
  optimal_shift=kron(expm(-1i*(pi/4)*creation(op_cutoff)*creation(op_cutoff)'),identity(op_cutoff));

  % Effect of the 50:50 beam splitter
  proj_columns=optimal_shift'*beam_splitter(op_cutoff)*proj_columns;

elseif pom_choice==4
% 4) Undoing the preparation of the initial state + photon counting

  if state_choice==1
  else
    error('This POM is only available for coherent states.')
  end

  % Observable quantity (number of photons at each port)
  observable=kron(creation(op_cutoff)*creation(op_cutoff)',creation(op_cutoff)*creation(op_cutoff)');
  [proj_columns,outcomes_temp]=eig(full(observable));
  outcomes=zeros(1,length(outcomes_temp));
  for x=1:length(outcomes_temp)
    outcomes(x)=outcomes_temp(x,x);
  end
  
  % Extra phase shifts
  optimal_shift=sparse(expm(-1i*pi*j3schwinger(op_cutoff)));

  % Unitary transformations to undo the preparation of the state
  bs=sparse(beam_splitter(op_cutoff)');
  cs_undo=sparse(kron(displacement(op_cutoff,sqrt(2)),identity(op_cutoff)));
  combined=cs_undo*bs*optimal_shift;
  proj_columns=combined'*proj_columns;

elseif pom_choice==5
  % 5) 50:50 beam splitter + parity measurements
  
  % Observable quantity (parity of the number of photons at each port)
  paritya=sparse(kron(identity(op_cutoff),(-1)^(full(creation(op_cutoff)*creation(op_cutoff)'))));
  parityb=sparse(kron((-1)^(full(creation(op_cutoff)*creation(op_cutoff)')),identity(op_cutoff)));
  observable=full(paritya*parityb);
  [proj_columns,outcomes_temp]=eig(full(observable));
  outcomes=zeros(1,length(outcomes_temp));
  for x=1:length(outcomes_temp)
    outcomes(x)=outcomes_temp(x,x);
  end

  % Extra phase shift
  odd_shift=kron(identity(op_cutoff),expm(1i*(pi/2)*creation(op_cutoff)*creation(op_cutoff)'));
  even_shift=kron(identity(op_cutoff),expm(1i*(pi/4)*creation(op_cutoff)*creation(op_cutoff)'));
  if state_choice==1
    optimal_shift=odd_shift;
  else
  optimal_shift=even_shift;
  end
  
  % Effect of the 50:50 beam splitter
  proj_columns=optimal_shift'*beam_splitter(op_cutoff)*proj_columns;

end
end
