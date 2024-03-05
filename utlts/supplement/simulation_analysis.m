function [psegs, Lams, ages, sub_ids] = save_outputs_analysis(data_dir,subs,correct_T,ages,Neff,norm_thresh)


TOTAL_SUBS = length(subs);
sub_ids = strings(1,TOTAL_SUBS); 
sub_ages = NaN(1,TOTAL_SUBS);    

Lams = NaN(1,TOTAL_SUBS); %do NaN mean to get rid of NaNs   
psegs = NaN(1,TOTAL_SUBS);    


[f2, norm, ~] = Ising_model(Neff);


for s=1:TOTAL_SUBS
	sub = subs{s};

	filename = strcat(data_dir,'/sim-',sub,'.csv'); 
    
	if isfile(filename)
		raw_data = readmatrix(filename);

		sz=size(raw_data);
		T=sz(1);
		num_regions=sz(2);

		if T == correct_T 
			sub_ids(s) = sub;
			sub_ages(s) = ages(s); 
			
			[Lams(s), psegs(s)] = fit_model(raw_data, correct_T, num_regions, Neff, norm_thresh, f2, norm);

		end

	end
end
end

function [Lam, pseg] = fit_model(raw_data, T, num_regions, Neff, norm_thresh, f2, norm)
	%synch=sum(Isingify2(T,raw_data),2)/num_regions;
	synch=sum(raw_data,2)/num_regions;	%only difference with fMRI data (no need to binarize)
		
	m2=mean(synch.^2);
	lam=fzero(@(lambda) f2(double(lambda))/norm(double(lambda)) - m2,  .00001, optimset('TolX',1*10^-10));                       %ML estimation
	Lamc=1/(2*Neff);
	Lam=(lam-Lamc)/Lamc;

%	m2_fit = f2(double(lam))/norm(double(lam))	#equals to m2

	pseg=sum(abs(synch)< norm_thresh )/T;
end
