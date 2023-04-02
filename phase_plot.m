for r1 = [3 5 10]
load(sprintf('phase_n300_c_r%i.mat',r1));
load(sprintf('phase_n300_f_r%i.mat',r1));

types = {'rc','fc','rf','ff'};
for t = 1:4
    type = types{t};
    if strcmp(type,'rc')
        errors = rc_error;
    elseif strcmp(type,'fc')
        errors = fc_error;
    elseif strcmp(type,'rf')
        errors = rf_error;
    else
        errors = ff_error;
    end
success = abs(errors) < 0.0001;
phase = flipud(mean(success,3));
figure
colormap('gray')
imagesc(phase)
if r1 == 10
    set(gca, 'XTick', 1:1:5, 'XTickLabel', 0.2:0.2:1)
else
    set(gca, 'XTick', 2:2:10, 'XTickLabel', 0.4:0.4:2)
end
set(gca, 'YTick', 1:5:31, 'YTickLabel', 0.6:-0.1:0)
set(gca, 'FontSize', 18)
xlabel('Sampling Constant', 'Interpreter','latex','Fontsize',24)
ylabel('Corruption Rate $(\%)$', 'Interpreter','latex','Fontsize',24)

fname_out = sprintf('results/phase_new_%s_r%d',type,r1);
saveas(gcf,fname_out,'eps')
end
end

