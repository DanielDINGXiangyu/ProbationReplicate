function [] = plot_Para(a,RRA,IES,name)
%This function compare different parameter setting

    aa = picker([a]);
    legend_str = {}; 
    colDif = 1/length(IES);
    
    figure(1)
    set(gcf,'position',[0,0,1500,500])
    subplot(1,3,1)
    for ii=1:length(IES)
        legend_str = [legend_str; sprintf('γ=%d ψ=%.1f ', RRA(ii), IES(ii))];
        load(sprintf('Solution/LowBound/nocorr_rra%02dies%02d.mat', RRA(ii), round(IES(ii)*10)));
        plot(B./g(1, aa), bp(:, aa), 'Color', [colDif*ii 1-colDif*ii 1],'Linewidth',1.5)
        hold on
    end
    hold off
    ylim([-1.1,-0.4]) 
    xlim([-1.1,-0.4]) 
    legend(legend_str,'Location','southeast')
    xlabel('Bond/Current Trend: $\frac{b_t}{\Gamma_{t}}$','interpreter', 'latex')
    ylabel("DE Bond Policy: $\frac{b_{t+1}}{\Gamma_{t}}$",'interpreter', 'latex')
    title("DE Bond Policy")
   
    
    %-------------------
    subplot(1,3,2)
    legend_str = {}; 
    for ii=1:length(IES)
        legend_str = [legend_str; sprintf('γ=%d ψ=%.1f ', RRA(ii), IES(ii))];
        load(sprintf('Solution/LowBound/nocorr_rra%02dies%02d.mat', RRA(ii), round(IES(ii)*10)));
        plot(B./g(1, aa), bpSP(:, aa), 'Color', [colDif*ii 1-colDif*ii 1],'Linewidth',1.5)
        hold on
    end
    hold off
    ylim([-1.1,-0.4]) 
    xlim([-1.1,-0.4]) 
    legend(legend_str,'Location','southeast')
    xlabel('Bond/Current Trend: $\frac{b_t}{\Gamma_{t}}$','interpreter', 'latex')
    ylabel("SP Bond Policy: $\frac{b_{t+1}}{\Gamma_{t}}$",'interpreter', 'latex')
    title("SP Bond Policy")
    
    
    %----------------------
    subplot(1,3,3)
    legend_str = {}; 
    for ii=1:length(IES)
        legend_str = [legend_str; sprintf('γ=%d ψ=%.1f ', RRA(ii), IES(ii))];
        load(sprintf('Solution/LowBound/nocorr_rra%02dies%02d.mat', RRA(ii), round(IES(ii)*10)));
        plot(B./g(1, aa), Tau(:, aa).*100, 'Color', [colDif*ii 1-colDif*ii  1],'Linewidth',1.5)
        hold on
    end
    hold off
    ylim([0,20]) 
    xlim([-1.1,-0.4])  
    legend(legend_str,'Location','southeast')
    xlabel('Bond/Current Trend: $\frac{b_t}{\Gamma_{t}}$','interpreter', 'latex')
    ylabel("Macroprudential Tax $\tau$ (\%)",'interpreter', 'latex')
    title("Optimal Tax")
    
    saveas(gcf,append(append('Plot/1tax/',name),'.png'))

end

