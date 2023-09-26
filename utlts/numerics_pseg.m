function [psegs]=numerics_pseg(Neff,seg_thresh,Lamage)
%% Pseg numerically calculated from Lambda

k_seg=(Neff-seg_thresh)/2:1:(Neff+seg_thresh)/2;
v_seg=(2*k_seg-Neff)/Neff;
vv_seg=v_seg.^2;
nck_seg=zeros(1,seg_thresh+1);
warning('off','all'); %choosek gives a warning
for n=1:seg_thresh+1
    nck_seg(n)=nchoosek(Neff,k_seg(n));
end

k=0:1:Neff;
v=(2*k-Neff)/Neff;
vv=v.^2;
nck=zeros(1,Neff+1);
warning('off','all'); %choosek gives a warning
for n=1:Neff+1
    nck(n)=nchoosek(Neff,k(n));
end


TOTAL_SUBS = length(Lamage);
psegs = ones(1,TOTAL_SUBS);
for i=1:TOTAL_SUBS
    
    Lamc=1/(2*Neff);
    lambda = (Lamage(i)+1)*Lamc;    %not sure why this gives different result than if reexpress everything in terms of Lambda, issue with 0s?

    amt_seg = sum(nck_seg.*exp(lambda.*vv_seg*Neff^2));
    amt_total = sum(nck.*exp(lambda.*vv*Neff^2));
    
    psegs(i)=amt_seg/amt_total; 
end
