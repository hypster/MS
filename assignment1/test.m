plot(dat0);
saveas(gcf, 'test.png');
histfit(dat0, [],'exponential');
saveas(gcf, 'histogram_exponential.png');
histfit(dat0, [],'gamma');
saveas(gcf,'histogram_gamma.png');




% hist(dat0);