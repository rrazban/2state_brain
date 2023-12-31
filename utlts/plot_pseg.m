function [total, hl, rho, pval]= plot_pseg(Edges, pseg, Sub_Ages, color, leg_lab)

    inds = ~isnan(pseg);	%remove Nans which correspond to unphysical value of error in reading in file
    pseg=pseg(inds);
    Sub_Ages=Sub_Ages(inds);

    N = length(Edges)-1;
    Q25e=zeros(1,N);
    Q75e=zeros(1,N);
    Qmed=zeros(1,N);

    total = 0;
    bin_size=zeros(1,N);
    for i=1:N

    	q=find(Sub_Ages>=Edges(i)&Sub_Ages<Edges(i+1));          %calculates the medians and quartiles (error bars) of each age group 
        sss=pseg(q);
        bin_size(i) = length(sss);
        total=total+length(sss);
        Qmed(i)=median(sss);
            
        stderror = std(sss)/sqrt(length(sss));
        Q25e(i)=-stderror;
        Q75e(i)=stderror;

	%quartiles intead
	%Q25e(i)=Qmed(i)-quantile(sss,.25);
	%Q75e(i)=quantile(sss,.75)-Qmed(i);
     end
     
     [rho,pval]=corr(Sub_Ages,pseg','Type','Spearman');
     str = strcat('\rho=', num2str(rho,'%.2f'),' (', num2str(pval,2), ') \newlineN=', string(total));	%some reason doesnt work if put direction after 'DisplayName'

     hl=errorbar(Edges(1:N)+2.5,Qmed,Q25e,Q75e,'DisplayName',str);

     hl.Marker='o';
     hl.MarkerEdgeColor= [.2 .2 .2];
     hl.MarkerFaceColor= color;
     hl.MarkerSize=10;  
     hl.LineWidth=1;
     hl.Color=color;
     
     
     hAxis=gca;
     hAxis.TickLength=[.04 .04];
     hAxis.XMinorTick='off';
     hAxis.YMinorTick='off';
     hAxis.LineWidth=1;
     hAxis.FontSize = 18;%default
     set(gca,'box','off');   %remove axes on top and right
     hold on
end
