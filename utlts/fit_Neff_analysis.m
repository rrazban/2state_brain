function [sub_ages, Neffs, m2s, lams] = fit_Neff_analysis(data_dir,subs,correct_T, ages)

TOTAL_SUBS = length(subs);
m4_accuracy = NaN(1,TOTAL_SUBS);    
m2s = NaN(1,TOTAL_SUBS);    
Neffs = NaN(1,TOTAL_SUBS);    
lams = NaN(1,TOTAL_SUBS);    
sub_ages = NaN(1,TOTAL_SUBS);    

possible_Neffs = 4:1:500;	%4 is the lowest you can go w/o errors being thrown, but lams are usually negative! start search at Neff~10 so that lam>0
%possible_Neffs = 4:1:100;	%runs faster but coverage is poorer since only search up to Neff=100	
min_n = 1;
max_n = length(possible_Neffs);

delete(gcp('nocreate'));	%delete any preexisting parallel session
parpool(3);	%3 parallel workers

%for s=1:TOTAL_SUBS
parfor s=1:TOTAL_SUBS
	sub = subs{s}%;
	
	filename = get_filename(data_dir, sub);

	if isfile(filename)
		if contains(filename, '.mat')	
			raw_data = load(filename).raw_data;
		else
			raw_data = readmatrix(filename);
		end	

		sz=size(raw_data);
		T=sz(1);
		num_regions=sz(2);	%also known as N

		if T == correct_T 
			sub_ages(s) = ages(s); 

			poss_m4_accs = NaN(1, length(possible_Neffs));
			poss_lams = NaN(1, length(possible_Neffs));
			poss_m2s = NaN(1, length(possible_Neffs));

			for n=min_n:max_n
				Neff = possible_Neffs(n);
				[f2, norm, f4] = Ising_model(Neff);
				[poss_m4_accs(n), poss_lams(n), poss_m2s(n)] = fit_model(raw_data, T, num_regions, f2, norm, f4, Neff);
			end
			[M,I]=min(abs(poss_m4_accs));	%can be negative and positive, thats why take abs!
			if I==min_n || I==max_n
            			disp(strcat(sub,' is at the Neff edge! value not recorded: ', num2str(possible_Neffs(I))))
				m2s(s)=poss_m2s(I);
			else
				m4_accuracy(s)=poss_m4_accs(I);
				lams(s)=poss_lams(I);
				Neffs(s)=possible_Neffs(I);

				m2s(s)=poss_m2s(I);
			end
		end
	end
end
end


function [m4_accuracy, lam, m2] = fit_model(raw_data, T, num_regions, f2, norm, f4, Neff)

	synch=sum(Isingify2(T,raw_data),2)/num_regions;
	m2=mean(synch.^2);
	
	lam=fzero(@(lambda) f2(double(lambda))/norm(double(lambda)) - m2,.00001, optimset('TolX',1*10^-10));                       %ML estimation

	m4_data=mean(synch.^4);
	m4_fit = f4(double(lam))/norm(double(lam));
	m4_accuracy = m4_fit-m4_data;

	%dont move this! cuz m4_fit depends on correct lam scaling
	lam=lam*Neff^2; 
end
