preNest(:,:,1) = fixShortNanGaps(preNest(:,:,1),10);
preNest(:,:,2) = fixShortNanGaps(preNest(:,:,2),10);

%%
theta = preNest(:,:,3);
rho = ones(size(theta));
[x y] = pol2cart(theta, rho);
scFact = 0.01;
x = x.*scFact;
y = y.*scFact;
xl = [-0.01 0.23];
yl = [-0.01 0.2];

for i = 2000:3000
    plot(preNest(i,:,1), preNest(i,:,2), 'ro');
    hold on
    plotbroodGreyTrans(broodPre, 200, 0.5, 0.2);
    plot(preNest(i,:,1), preNest(i,:,2), 'ro');
    
    quiver(preNest(i,:,1), preNest(i,:,2), x(i,:),y(i,:), 0);
    hold off
    xlim(xl);
    ylim(yl);
    drawnow
end
