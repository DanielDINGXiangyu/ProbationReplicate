function [] = plot_VolTax_vertical(a,b,c,d, RRA,IES,name)
% Plot figure: volatility tax across z and g volatility regime

    aa = picker([a]); % low Z vol
    bb = picker([b]); % high Z vol
    cc = picker([c]); % low G vol 
    dd = picker([d]); % high G vol
    %-------------------------------------
    set(gcf,'position',[0,0,length(IES)*200,1200])
   
    
    for ii=1:length(IES)
        %% aa and bb are compare Low and High Z vol
        subplot(4,ceil(length(IES)/2),ii)
        load(sprintf('Solution/LowBound/nocorr_rra%02dies%02d.mat', RRA(ii), round(IES(ii)*10)));
        plot(B./g(1, aa), Tau(:, aa).*100, 'r','Linewidth',1.5)
        hold on
        plot(B./g(1, bb), Tau(:, bb).*100, 'b','Linewidth',1.5)
        hold off
        ylim([0,20]) 
        xlim([-1.1,-0.4]) 
        legend({'Low Z Vol: SP', 'High Z Vol: SP'},'Location','southeast')
        xlabel('Bond/Current Trend: $\frac{b_t}{\Gamma_{t}}$','interpreter', 'latex')
        ylabel("Macroprudential Tax $\tau$ (\%) ",'interpreter', 'latex')
        title(sprintf('Across Z Vol; γ=%d ψ=%.1f ', RRA(ii), IES(ii)))

        %% cc and dd are compare Low and High G vol
        subplot(4,ceil(length(IES)/2),ii+length(IES))
        plot(B./g(1, cc), Tau(:, cc).*100, 'r','Linewidth',1.5)
        hold on
        plot(B./g(1, dd), Tau(:, dd).*100, 'b','Linewidth',1.5)
        hold off
        ylim([0,20]) 
        xlim([-1.1,-0.4]) 
        legend({'Low G Vol: SP', 'High G Vol: SP'},'Location','southeast')
        xlabel('Bond/Current Trend: $\frac{b_t}{\Gamma_{t}}$','interpreter', 'latex')
        ylabel("Macroprudential Tax $\tau$ (\%) ",'interpreter', 'latex')
        title(sprintf('Across G Vol; γ=%d ψ=%.1f ', RRA(ii), IES(ii)))

    end

    saveas(gcf,append(append('Plot/4taxRobustness/',name),'_NoCorr.png'))
end
