function [m4_data, m4_fit, lams] = identify_Neff_analysis(data_dir,subs,correct_T,Neff)


TOTAL_SUBS = length(subs);
m4_data = NaN(1,TOTAL_SUBS); %do NaN mean to get rid of NaNs   
m4_fit = NaN(1,TOTAL_SUBS);    
m2s = NaN(1,TOTAL_SUBS);    
lams = NaN(1,TOTAL_SUBS);    


[f2, norm, f4] = Ising_model(Neff);

delete(gcp('nocreate'));	%delete any preexisting parallel session
parpool(3);	%3 parallel workers

%for s=1:TOTAL_SUBS	
parfor s=1:TOTAL_SUBS
	sub = subs{s};
	
	filename = get_filename(data_dir, sub);
   
	if isfile(filename)
		if contains(filename, '.mat')
			raw_data = load(filename).raw_data;
		else
			raw_data = readmatrix(filename);
		end	

		%for UKB, remove first line and first column cuz they are just labels
		if contains(data_dir, 'ukb')
			raw_data(1, :) = [];
			raw_data(:, 1) = [];
		end

		sz=size(raw_data);
		T=sz(1);
		num_regions=sz(2);	%also known as N

		if T == correct_T 
			[m4_data(s), m4_fit(s), lams(s), m2s(s)] = fit_model(raw_data, T, num_regions, f2, norm, f4, Neff);
		end
	end
end

end


function [m4_data, m4_fit, lam, m2] = fit_model(raw_data, T, num_regions, f2, norm, f4, Neff)
		synch=sum(Isingify2(T,raw_data),2)/num_regions;
		
		m2=mean(synch.^2);
		Lamc=1/(2*Neff);
		lam=fzero(@(lambda) f2(double(lambda))/norm(double(lambda)) - m2, 0.8*Lamc, optimset('TolX',1*10^-10));	%set initial guess at Lam=-0.2	#-0.1 sometimes gives errors at various Neff values

		m2_fit = f2(double(lam))/norm(double(lam));	%need double precission, or else overflows for hcp when Neff around >100, end up getting lots of NaNs

		m4_data=mean(synch.^4);
		m4_fit = f4(double(lam))/norm(double(lam));
end
