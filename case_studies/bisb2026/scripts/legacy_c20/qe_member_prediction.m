function [Y,B,C,J] = qe_member_prediction(E,q,p,n,model_name,ampunit)
% p=[r,B(1:Q), {Eleft,Eright,width,A(1:Q)} x n]. One regional scale.
E=E(:); q=q(:).'; Q=numel(q); ne=numel(E); beta=(q-min(q))/(max(q)-min(q));
shape=(E/1000).^(-p(1)); B=shape*p(2:Q+1); C=zeros(ne,Q,n);
wantJ=nargout>3; J=[];
if wantJ
 J=zeros(ne*Q,numel(p)); J(:,1)=reshape(-log(E/1000).*B,[],1);
 for k=1:Q, J((k-1)*ne+(1:ne),k+1)=shape; end
end
for j=1:n
 base=1+Q+(j-1)*(3+Q); center=1000*((1-beta)*p(base+1)+beta*p(base+2)); width=1000*p(base+3);
 for k=1:Q
  amplitude=p(base+3+k)*ampunit; d=E-center(k);
  if strcmp(model_name,'lorentz_symmetric')
   D=d.^2+width^2/4; unit=width/(2*pi)./D;
   val=amplitude*unit; dc=2*d./D.*val;
   dw=amplitude/(2*pi)*(d.^2-width^2/4)./D.^2;
  else
   D=(E.^2-center(k)^2).^2+E.^2*width^2; unit=E*width./D;
   val=amplitude*unit; dc=val.*(4*center(k)*(E.^2-center(k)^2))./D;
   dw=amplitude*E./D-val.*(2*E.^2*width)./D;
  end
  C(:,k,j)=val;
  if wantJ
   rows=(k-1)*ne+(1:ne);
   J(rows,base+1)=dc*1000*(1-beta(k)); J(rows,base+2)=dc*1000*beta(k);
   J(rows,base+3)=dw*1000; J(rows,base+3+k)=unit*ampunit;
  end
 end
end
Y=B+sum(C,3);
end
