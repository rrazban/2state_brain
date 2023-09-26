function [f2, normalization, f4] = Ising_model(Neff)
	k=0:1:Neff;
	v=(2*k-Neff)/Neff;
	vv=v.^2;
	vvvv=vv.^2;

	nck=zeros(1,Neff+1);
	warning('off','all'); %choosek gives a warning about inaccuracy about 10^15


	for n=1:Neff+1
		nck(n)=nchoosek(Neff,k(n));	%at least correct order of magnitude by taking 10^above 15 power
%		nck(n)=nchoosek(uint64(Neff),k(n));	%this arbitrarily truncates sum	%need uint64 so no overflow for large N
	end
%rescale \lambda by Neff^2 to follow convention from Weistuch et al. 2021
%fitting works fine either way

	f2=@(lambda) sum((vv).*nck.*exp(lambda.*vv*Neff^2));   %mean s^2 (variance of s)
%	f2=@(lambda) sum((vv).*nck.*exp(lambda.*vv));   %mean s^2 (variance of s)

	normalization=@(lambda) sum(nck.*exp(lambda.*vv*Neff^2));         %normalization         
%	normalization=@(lambda) sum(nck.*exp(lambda.*vv));         %normalization         

	f4=@(lambda) sum(vvvv.*nck.*exp(lambda.*vv*Neff^2)); 
%	f4=@(lambda) sum(vvvv.*nck.*exp(lambda.*vv)); 
end
