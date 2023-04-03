function alpha=friction_factor(mue,nn,free_sld,H,C)
% this function calculate the friction factor for each free sliding node

alpha=ones(nn,1);
% calculat angle 
% angle of circluar strings
phi=zeros(nn,1);
for i=1:numel(free_sld)
    h_1=H(:,find(C(:,free_sld(i))==1));
    h_2=H(:,find(C(:,free_sld(i))==-1));
phi(free_sld(i))=acos(-h_1'*h_2/norm(h_1)/norm(h_2));% approximation not accurate
end
alpha(free_sld)=exp(-mue*phi(free_sld));


end