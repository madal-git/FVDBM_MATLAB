figure;
for r=1:M
    P=CELL{r};
    plot(P{22},P{23},'black',P{24},P{25},'black',P{26},P{27},'black','linewidth',0.5);
    hold on
end


% for r=1:M
%     P=CELL{r};
%     plot(P{22},P{23},P{24},P{25},P{26},P{27},'black', 'linewidth',1);
%     hold on;
% %     plot(P{24},P{25},'black', 'linewidth',1);
% %     hold on;
% %     plot(P{26},P{27},'black', 'linewidth',1);
% %     hold on;
% end
axis equal tight