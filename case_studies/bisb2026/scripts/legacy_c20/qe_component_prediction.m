function [total,bg,components] = qe_component_prediction(E,p,model_name,n,ampunit)
% Solver coordinates; output scaled ordinate, never extra dE.
model=peak_models(model_name);
bg=p(1)*(E/1000).^(-p(2));
if numel(p)>2+3*n, bg=bg+p(end); end
components=zeros(numel(E),n);
for j=1:n
 ii=3+(j-1)*3;
 components(:,j)=model.model_fn(p(ii)*1000,p(ii+1)*1000,p(ii+2)*ampunit,E);
end
total=bg+sum(components,2);
end
