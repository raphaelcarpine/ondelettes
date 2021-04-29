figure;

%%
c = colororder;
c = [c; rand(2, 3)]

colororder(c);


%%
hold on
for r=1:9
    x = linspace(0,r,500);
    y = sqrt(r.^2-x.^2);
    plot(x,y,'LineWidth',15);
end