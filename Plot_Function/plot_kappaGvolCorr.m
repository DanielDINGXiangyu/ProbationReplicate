
function [] = plot_kappaGvolCorr(a,b,RRA,IES,name)
%% 1. Kappa correlate with local growth volatility

    aa = picker([a]);
    bb = picker([b]);

    figure(1)
    set(gcf,'position',[0,0,length(IES)*500,1000])


    for ii=1:length(IES)
        subplot(2,length(IES),ii)
        load(sprintf('Solution/LowBound/kappagcorr_rra%02dies%02d_spreadkappa_lowsticky.mat', RRA(ii), round(IES(ii)*10)));
        plot(B./g(1, aa), bp(:, aa), 'r','Linewidth',1.5)
        hold on
        plot(B./g(1, aa), bpSP(:, aa), 'r--','Linewidth',1.5)
        hold on
        plot(B./g(1, bb), bp(:, bb), 'b','Linewidth',1.5)
        hold on
        plot(B./g(1, bb), bpSP(:, bb), 'b--','Linewidth',1.5)
        hold off
        ylim([-1.1,-0.4]) 
        xlim([-1.1,-0.4]) 
        legend({'Low G Vol: DE', 'Low G Vol: SP', 'High G Vol: DE','High G Vol: SP'},'Location','southeast')
        xlabel('Bond/Current Trend: $\frac{b_t}{\Gamma_{t}}$','interpreter', 'latex')
        ylabel("Detrended Bond Policy: $\frac{b_{t+1}}{\Gamma_{t}}$",'interpreter', 'latex')
        title(sprintf('Across G Vol; γ=%d ψ=%.1f ', RRA(ii), IES(ii)))


        subplot(2,length(IES),ii+length(IES))
        plot(B./g(1, aa), Tau(:, aa).*100, 'r','Linewidth',1.5)
        hold on
        plot(B./g(1, bb), Tau(:, bb).*100, 'b','Linewidth',1.5)
        hold off
        ylim([0,20]) 
        xlim([-1.1,-0.4]) 
        legend({'Low G Vol: SP', 'High G Vol: SP'},'Location','southeast')
        xlabel('Bond/Current Trend: $\frac{b_t}{\Gamma_{t}}$','interpreter', 'latex')
        ylabel("Macroprudential Tax $\tau$ (\%) ",'interpreter', 'latex')
        title(sprintf('Across G Vol; γ=%d ψ=%.1f ', RRA(ii), IES(ii)))

    end

    saveas(gcf,append(append('Plot/2kappa/',name),'_withCorr.png'))



    %-------------------------------------
    %% 2. Compare to benchmark case without any correlation

    figure(2)
    set(gcf,'position',[0,0,length(IES)*500,1000])

    for ii=1:length(IES)
        subplot(2,length(IES),ii)
        load(sprintf('Solution/LowBound/nocorr_rra%02dies%02d.mat', RRA(ii), round(IES(ii)*10)));
        plot(B./g(1, aa), bp(:, aa), 'r','Linewidth',1.5)
        hold on
        plot(B./g(1, aa), bpSP(:, aa), 'r--','Linewidth',1.5)
        hold on
        plot(B./g(1, bb), bp(:, bb), 'b','Linewidth',1.5)
        hold on
        plot(B./g(1, bb), bpSP(:, bb), 'b--','Linewidth',1.5)
        hold off
        ylim([-1.1,-0.4]) 
        xlim([-1.1,-0.4]) 
        legend({'Low G Vol: DE', 'Low G Vol: SP', 'High G Vol: DE','High G Vol: SP'},'Location','southeast')
        xlabel('Bond/Current Trend: $\frac{b_t}{\Gamma_{t}}$','interpreter', 'latex')
        ylabel("Detrended Bond Policy: $\frac{b_{t+1}}{\Gamma_{t}}$",'interpreter', 'latex')
        title(sprintf('Across G Vol; γ=%d ψ=%.1f ', RRA(ii), IES(ii)))


        subplot(2,length(IES),ii+length(IES))
        plot(B./g(1, aa), Tau(:, aa).*100, 'r','Linewidth',1.5)
        hold on
        plot(B./g(1, bb), Tau(:, bb).*100, 'b','Linewidth',1.5)
        hold off
        ylim([0,20]) 
        xlim([-1.1,-0.4]) 
        legend({'Low G Vol: SP', 'High G Vol: SP'},'Location','southeast')
        xlabel('Bond/Current Trend: $\frac{b_t}{\Gamma_{t}}$','interpreter', 'latex')
        ylabel("Macroprudential Tax $\tau$ (\%) ",'interpreter', 'latex')
        title(sprintf('Across G Vol; γ=%d ψ=%.1f ', RRA(ii), IES(ii)))

    end

    saveas(gcf,append(append('Plot/2kappa/',name),'_NoCorr.png'))
end
