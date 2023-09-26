function [psegs, Lams, ages, sub_ids] = save_outputs_analysis(data_dir,subs,correct_T,ages,Neff,norm_thresh)


TOTAL_SUBS = length(subs);

[f2, norm, ~] = Ising_model(Neff);


for s=1:TOTAL_SUBS
	sub = subs{s}

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
			
			synch=sum(Isingify2(T,raw_data),2)/num_regions;
			histogram(synch,10,'BinLimits',[-1,1],'normalization','pdf','DisplayStyle','stairs','LineWidth',4)
			hold on

			[f2_std, norm_std, ~] = Ising_model(num_regions);
			lam_std = fit_model(synch, f2_std, norm_std);
			plot_p_s(lam_std, num_regions)
			
			
			lam = fit_model(synch, f2, norm);
			plot_p_s(lam, Neff)
			
			xlabel('synchrony $s$', 'Interpreter','latex');
			ylabel('$P(s)$','Interpreter','latex');
			hAxis=gca;
			hAxis.LineWidth=1;
			hAxis.FontSize = 18;
			legend('human subject (300 regions)', 'Ising model (300 regions)', 'Ising model ('+string(Neff)+' effective regions)','Location', 'southeast', 'FontSize',12)
			pause
			hold off
		end

	end
end
end

function lam = fit_model(synch, f2, norm)
		
	m2=mean(synch.^2);
	lam=fzero(@(lambda) f2(double(lambda))/norm(double(lambda)) - m2,  .00001, optimset('TolX',1*10^-10));                       %ML estimation
end

function plot_p_s(lam, Neff)
	k=0:1:Neff;
	v=(2*k-Neff)/Neff;
	vv=v.^2;
	vvvv=vv.^2;

	nck=zeros(1,Neff+1);
	warning('off','all'); %choosek gives a warning about inaccuracy about 10^15


	for n=1:Neff+1
		nck(n)=nchoosek(Neff,k(n));
	end



	pfit=nck.*exp(lam*vv*Neff^2);

	F = griddedInterpolant(v,pfit);
	fun = @(t) F(t);
	Z = integral(fun, v(1), v(end));    %normalization to be able to directly compare P(s) with histogram
	plot(v,pfit/Z,'LineWidth',4);  %plot wrt s^2, the variable whose mean we are fitting
	hold on
end
