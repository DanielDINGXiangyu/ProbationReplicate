function [] = plot_dePolicy(a,b,c,RRA,IES,name)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    aa = picker([a]);
    bb = picker([b]);
    cc = picker([c]);
    
    if abs(a-c)> 900
        str = "Z Level";
    elseif abs(a-c)> 90
        str = "Z Vol";
    elseif abs(a-c)>9
        str = "G Level";
    else
        str = "G Vol";
    end
        
    
    set(gcf,'position',[0,0,length(IES)*500,500])
    for ii=1:length(IES)
        subplot(1,length(IES),ii)
        load(sprintf('Solution/LowBound/nocorr_rra%02dies%02d.mat', RRA(ii), round(IES(ii)*10)));
        plot(B./g(1, aa), bp(:, aa), 'r','Linewidth',1.5)
        hold on
        plot(B./g(1, bb), bp(:, bb), 'g','Linewidth',1.5)
        hold on
        plot(B./g(1, cc), bp(:, cc), 'b','Linewidth',1.5)
        hold off
        ylim([-1.1,-0.4]) 
        xlim([-1.1,-0.4]) 
        legend({sprintf('Low %s: DE', str), sprintf('Mid %s: DE', str), sprintf('High %s: DE', str)},'Location','southeast')
        xlabel('Bond/Current Trend: $\frac{b_t}{\Gamma_{t}}$','interpreter', 'latex')
        ylabel("Detrended Bond Policy: $\frac{b_{t+1}}{\Gamma_{t}}$",'interpreter', 'latex')
        title(sprintf('Across %s; γ=%d ψ=%.1f ',str, RRA(ii), IES(ii)))  
    end
    
    saveas(gcf,append(append('Plot/0ShockParameter/',name),'_DE.png'))
    
end

