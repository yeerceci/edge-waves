% figure
hold on
set(0, 'DefaultLineLineWidth', 1);
% l_plot( find(l_plot>0) ) = NaN;

if iBC == 0,
    plot(real(l_plot), w_plot, 'o')
%     xlim([-200  0])
    ylim([0 10])
%     legend('Real bathymetry')
elseif iBC == 0.1, % green cross
    plot(real(l_plot), w_plot, 'rx','color','[0.6350    0.0780    0.1840]')
elseif iBC == 0.4, % deep red cross
    plot(real(l_plot), w_plot, 'x', 'color','[0.4660    0.6740    0.1880]')
elseif iBC == 1,  % Balls - purple dashed
    plot(l_plot, w_plot,'--','color','[0.4940    0.1840    0.5560]')
elseif iBC == 4,  % wedge - blue solid
    plot(l_plot, w_plot,'color','[ 0    0.4470    0.7410]')
end
    
% legend('Wedge - analytical', 'Wedge - discrete','Basin - analytical', 'Basin - discrete')
% legend('Wedge - analytical', 'Wedge - discrete','Basin - analytical', 'Basin - discrete', 'Real bathymetry')
% legend('Wedge', 'Basin', 'Real bathymetry')
% xlim([-50  10])
% ylim([0  1.5])
xlabel('gl/f^2')
ylabel('\omega/f')
title('n = 0, l < 0 mode')
set(gca,'FontSize',14)
set(gcf,'Position',[0 0 400 600])
grid on