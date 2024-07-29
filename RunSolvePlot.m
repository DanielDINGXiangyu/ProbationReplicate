
%% Solve the model under different parameter 
solve_nocorr(8,1.5,'rra08ies15')
solve_nocorr(8,0.5,'rra08ies05')
solve_nocorr(2,1.5,'rra02ies15')
solve_nocorr(2,0.5,'rra02ies05')


%% Replicate Figure 1
plot_VolTax_vertical(3132,3332,3231,3233,[2 8 2 8], [0.5 0.5 1.5 1.5],sprintf('FigureVol_cycle%01dtrend%01d.mat',3,3));


%% Replicate Figure 2&3

plot_mechanism_numden(3132,3332, [2 8], [0.5 0.5],true,'zVol_freeaxis')
plot_mechanism_numden(3231,3233, [2 8], [0.5 0.5],true,'gVol_freeaxis')