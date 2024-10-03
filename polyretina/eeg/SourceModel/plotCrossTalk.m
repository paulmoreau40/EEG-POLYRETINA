function plotCrossTalk(CT, fields)
figure
hold on
imagesc(1:numel(fields),1:numel(fields),CT, [0,1]);
colormap('gray');
xline(3.5,'--r', 'LineWidth', 2);
yline(3.5,'--r', 'LineWidth', 2);
axis image ij
xticklabels(fields)
xlabel('Receiving area')
yticklabels(fields)
ylabel('Seed area')
colorbar;
end