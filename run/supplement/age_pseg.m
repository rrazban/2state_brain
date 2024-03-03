function age_pseg(dataset, parcellation, which)

[title_text, Edges, data_dir, ~, ~, ~, ~,~,~] = feval(dataset, parcellation); 

filename = strcat('output_',parcellation ,'.mat');

load(strcat(dataset,'/',filename),'subs','ages','psegs','Lams','pred_psegs')

inds=Lams<-1;
sum(inds)
ages(inds)=NaN;%NaNs get removed in plot_pseg
psegs(inds)=NaN;
Lams(inds)=NaN;


%%code to only include individuals of a certain sex
%info_file = 'CAMCAN_participant_data.tsv';
%info = tdfread(info_file);
%all_subs = cellstr(info.Observations);
%sex = info.gender_code;
%sex_subs = all_subs(sex==1);	%1 (male) or 2 (female) for CamCAN
%[C, ia, ib] = intersect(subs, sex_subs);
%ages=ages(ia);
%psegs=psegs(ia);
%Lams=Lams(ia);



if contains(which, 'lam')
	ys=Lams;
	y_text = '\Lambda';
	location='southwest';
else
	ys=psegs;
	y_text = '{\it P}_{seg}';
	location='southeast';
end

hl = scatter(ages, ys, 'filled');
hold on

[n_exp, data, rho, pval] = plot_pseg(Edges,ys,ages,'m','');%[0.3010 0.7450 0.9330], [0.8500 0.3250 0.0980]

inds = ~isnan(ys);
ys=ys(inds);
ages=ages(inds);


xlim([min(Edges)-1 max(Edges)+1])
title(title_text)
ylabel(y_text);
xlabel('age in years')

lgnd = legend([hl data],{strcat('\rho=', num2str(rho,'%.2f'),' (', num2str(pval,2), ') \newlineN=', string(n_exp)), 'binned data'}, 'Location', location)


end

