function identify_Neff(dataset, parcellation)
%%identify Neff hyperparameter based on best match with synchrony fourth moment
%%also can show <s^4> vs Lambda plot

%dataset: 'camcan', 'hcp', 'ukb'
%parcellation: '300', 'gm_voxel'

hold on

%possible_Neffs=[40];
%possible_Neffs = [20 25 30 35 40 45 50];	
possible_Neffs = [115 120 125 130 135];	

max_n = length(possible_Neffs);

rmses = ones(1,max_n);    

for n=1:max_n 
    tic

    Neff = possible_Neffs(n); 

    [title_text, ~, data_dir, ~, subs, correct_T, ~,~,~] = feval(dataset, parcellation);	%feval evaluates the function of the corresponding dataset to get general info of dataset

    [m4_data, m4_fit, lams] = identify_Neff_analysis(data_dir, subs, correct_T, Neff);
    toc     

    inds = lams<0;	%get rid of unphysical fits
    disp(['number of unphysical lambda fits: ', num2str(sum(inds))])
    
    lams(inds)=NaN;
    m4_data(inds)=NaN;
    m4_fit(inds)=NaN;

%    plot_Lam_s4(lams, m4_data, m4_fit, title_text, Neff)	%check out fourth moment and corresponding fit
 %   pause

    rmses(n) = sqrt(mean((m4_fit(:)-m4_data(:)).^2, "omitnan"))
end
possible_Neffs

xlabel = 'effective number of regions{\it N}_{eff}';
ylabel = '\langle s^4\rangle RMSE';

plot_Neff(possible_Neffs, rmses, xlabel, ylabel, title_text, nnz(~isnan(m4_fit)))
end


function plot_Lam_s4(lams, m4_data, m4_fit, title_text, Neff)
%extra function for plotting <s^4> vs Lambda, make sure to uncomment within for loop above

	Lamc=1/(2*Neff);
	Lams=(lams-Lamc)/Lamc;

	exp_data = scatter(Lams, m4_data, 100);


	[out,idx] = sort(Lams);
%	plot(0,0)	%get right color for Ising model with Neff in paper figure
	theory = plot(out, m4_fit(idx), 'LineWidth', 4);	%sort so that it is a nice smooth line

	xlabel('rescaled connection strength \Lambda')
	ylabel('\langle s^4 \rangle')

	title(title_text)

	legend([exp_data theory], {'data (300 regions)', 'Ising model (' + string(Neff)+' effective regions)'},'Location', 'northwest', 'FontSize',14)

	hAxis=gca;
	hAxis.LineWidth=1;
	hAxis.FontSize = 18;
end
