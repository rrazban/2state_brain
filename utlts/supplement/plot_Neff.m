function plot(x,y, xtext, ytext, title_text, num_subs)
	scatter(x,y, 100,'MarkerEdgeColor',[0.3010 0.7450 0.9330],'MarkerFaceColor',[0.3010 0.7450 0.9330])

	xlabel(xtext)
	ylabel(ytext)

%	xlim([min(x)-2 max(x)+2])	
	[rho,pval] = corr(x',y', 'Type', 'Spearman', 'rows', 'complete')
	legend(strcat('data (N=',string(num_subs) ,')'), 'Location', 'northwest')

	hAxis=gca;
	title(hAxis, title_text, ' ')	%make sure title has room for yaxis overflow 
	hAxis.LineWidth=1;
	hAxis.FontSize = 18;
	
	ind = isnan(x);
	x(ind) = [];
	y(ind) = [];
	P = polyfit(x,y,1)
end
