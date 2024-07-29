function [] = plot_mechanism_numden(a,b, RRA,IES,freeaxis,name)
% Plot figure: world interest rate and growth volaitlity correlation

    aa = picker([a]); % low Z vol
    bb = picker([b]); % high Z vol
    %-------------------------------------
    set(gcf,'position',[0,0,length(IES)*500,1000])
    
    
    for ii=1:length(IES)
        %% plot numenator
        subplot(2,length(IES),ii)
        load(sprintf('Solution/LowBound/nocorr_rra%02dies%02d.mat', RRA(ii), round(IES(ii)*10)));
        plot(B./g(1, aa), Tau_num(:, aa), 'r','Linewidth',1.5)
        hold on
        plot(B./g(1, bb), Tau_num(:, bb), 'b','Linewidth',1.5)
        hold off
        if freeaxis==false
            ylim([0,0.36]) 
        end
        xlim([-1.1,-0.4]) 
        legend({'Low Vol: SP', 'High Vol: SP'},'Location','southeast')
        xlabel('Bond/Current Trend: $\frac{b_t}{\Gamma_{t}}$','interpreter', 'latex')
        ylabel("Tax Numerator ",'interpreter', 'latex')
        title(sprintf('Tax Numerator; γ=%d ψ=%.1f ', RRA(ii), IES(ii)))

        %% plot denomenator
        subplot(2,length(IES),ii+length(IES))
        plot(B./g(1, aa), Tau_den(:, aa), 'r','Linewidth',1.5)
        hold on
        plot(B./g(1, bb), Tau_den(:, bb), 'b','Linewidth',1.5)
        hold off
        if freeaxis==false
            ylim([0,5]) 
        end
        xlim([-1.1,-0.4]) 
        legend({'Low Vol: SP', 'High Vol: SP'},'Location','southeast')
        xlabel('Bond/Current Trend: $\frac{b_t}{\Gamma_{t}}$','interpreter', 'latex')
        ylabel("Tax Denominator ",'interpreter', 'latex')
        title(sprintf('Tax Denominator; γ=%d ψ=%.1f ', RRA(ii), IES(ii)))

    end

    saveas(gcf,append(append('Plot/5Mechanism/',name),'_NoCorr.png'))
end