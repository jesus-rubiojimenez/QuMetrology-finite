function [creat] = creation(dimension)
% Matrix representation of the creation operator, where ’dimension’ is
% the cutoff of the space.
creat=zeros(dimension,dimension);
for aa=1:dimension
  for bb=1:dimension
    if aa == (bb+1); creat(aa,bb)=sqrt(bb);
    else; creat(aa,bb)=0;
    end
  end
end
creat=sparse(creat); end
