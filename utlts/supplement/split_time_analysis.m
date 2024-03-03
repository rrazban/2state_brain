function [psegs, Lams, ages, sub_ids] = save_outputs_analysis(data_dir,subs,correct_T,ages,Neff,norm_thresh)


TOTAL_SUBS = length(subs);
sub_ids = strings(1,TOTAL_SUBS); 
sub_ages = NaN(1,TOTAL_SUBS);    

num_splits = 5;
Lams = NaN(num_splits,TOTAL_SUBS); %do NaN mean to get rid of NaNs   
psegs = NaN(num_splits,TOTAL_SUBS);    


[f2, norm, ~] = Ising_model(Neff);

%delete(gcp('nocreate'));	%delete any preexisting parallel session
%parpool(3);	%3 parallel workers


for s=1:TOTAL_SUBS
%parfor s=1:TOTAL_SUBS	%run code in parallel	%does not like the num_split for loop
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
		num_regions=sz(2);

		if T == correct_T 
			sub_ids(s) = sub;
			sub_ages(s) = ages(s); 

			t_split = fix(correct_T/num_splits); %fix forces division to be lower integer
			for i = 1:num_splits				
				data_subset = raw_data((i-1)*t_split+1:i*t_split, :);	%matlab splicing is inclusive on both sides
				[Lams(i,s), psegs(i,s)] = fit_model(data_subset, t_split, num_regions, Neff, norm_thresh, f2, norm);
			end
		end

	end
end
end

function [Lam, pseg] = fit_model(raw_data, T, num_regions, Neff, norm_thresh, f2, norm)
	synch=sum(Isingify2(T,raw_data),2)/num_regions;
		
	m2=mean(synch.^2);
	lam=fzero(@(lambda) f2(double(lambda))/norm(double(lambda)) - m2,  .00001, optimset('TolX',1*10^-10));                       %ML estimation
	Lamc=1/(2*Neff);
	Lam=(lam-Lamc)/Lamc;

%	m2_fit = f2(double(lam))/norm(double(lam))	#equals to m2

	pseg=sum(abs(synch)< norm_thresh )/T;
end
