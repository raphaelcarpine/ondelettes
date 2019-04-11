function niquistFrac()

clear all;
close all;


mu = 0.01;
omega0 = 1;
zeta0 = 0.0;
omega1 = omega0/(1+mu);
alpha = 1;
nZeta1 = 100;
Zeta10 = 0;
Zeta11 = 0.3;
Zeta1 = linspace(Zeta10, Zeta11, nZeta1);



%%
racines = nan(4, length(Zeta1));
for k = 1:length(Zeta1)
    zeta1 = Zeta1(k);
    racines(:,k) = polesSystFrac(mu, omega0, omega1, zeta0, zeta1, alpha);
end
racines = racines(~isnan(racines));
racines = racines(imag(racines) >= 0); %on garde seulement les poles à partie imaginaire positive pour l'affichage


%%
f = figure;
f.Position(2) = f.Position(2) - 0.5*f.Position(4);
f.Position(3:4) = f.Position(3:4) + 0.5*f.Position(3:4);


zeta1 = Zeta1(1);

%poles
ax = axes('Parent', f, 'Position', [0.05 0.25 0.9 0.7]);
hold(ax, 'on');
colors = get(ax, 'ColorOrder');
index = get(ax,'ColorOrderIndex');
poles = plot(real(racines)', imag(racines)', '.', 'Parent', ax, 'Color', colors(index, :));
xlabel('Re');
ylabel('Im');

racines = polesSystFrac(mu, omega0, omega1, zeta0, zeta1, alpha);
racines = racines(imag(racines) >= 0);
points = plot(real(racines), imag(racines), 'o', 'Parent', ax);
hold(ax, 'off');

grid(ax, 'on');



%%
b = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.05 0.14 0.45 0.03],'Style','slider',...
    'value',zeta1, 'min', Zeta10, 'max', Zeta11);
bgcolor = f.Color;
b2 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.31 0.09 0.15 0.05],'Style','edit',...
    'String',num2str(zeta1), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.2 0.09 0.1 0.05],'Style','text',...
    'String','zeta1', 'BackgroundColor',bgcolor);

b.Callback = @(es,ed) updateZeta(es.Value);
b2.Callback = @(es,ed) updateZeta(es.String);


bmu = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.01 0.01 0.09 0.04],'Style','edit',...
    'String',num2str(mu), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.01 0.05 0.09 0.04],'Style','text',...
    'String','mu', 'BackgroundColor',bgcolor);

bw = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.11 0.01 0.09 0.04],'Style','edit',...
    'String',num2str(omega1/omega0), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.11 0.05 0.09 0.04],'Style','text',...
    'String','w1/w0', 'BackgroundColor',bgcolor);

bzeta0 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.21 0.01 0.09 0.04],'Style','edit',...
    'String',num2str(zeta0), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.21 0.05 0.09 0.04],'Style','text',...
    'String','zeta0', 'BackgroundColor',bgcolor);

bZeta10 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.31 0.01 0.09 0.04],'Style','edit',...
    'String',num2str(Zeta10), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.31 0.05 0.09 0.04],'Style','text',...
    'String','Zeta1(1)', 'BackgroundColor',bgcolor);

bZeta11 = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.41 0.01 0.09 0.04],'Style','edit',...
    'String',num2str(Zeta11), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.41 0.05 0.09 0.04],'Style','text',...
    'String','Zeta1(end)', 'BackgroundColor',bgcolor);

balpha = uicontrol('Parent',f, 'Units', 'normalized','Position', [0.51 0.01 0.09 0.04],'Style','edit',...
    'String',num2str(alpha), 'BackgroundColor',bgcolor);
uicontrol('Parent',f, 'Units', 'normalized','Position', [0.51 0.05 0.09 0.04],'Style','text',...
    'String','alpha', 'BackgroundColor',bgcolor);

bmu.Callback = @(~,~) updateAll();
bw.Callback = @(~,~) updateAll();
bzeta0.Callback = @(~,~) updateAll();
bZeta10.Callback = @(~,~) updateAll();
bZeta11.Callback = @(~,~) updateAll();
balpha.Callback = @(~,~) updateAll();


%%

    function updateZeta(zeta)
        if nargin == 1
            if isnumeric(zeta)
                zeta1 = zeta;
                set(b2, 'String', num2str(zeta1));
            else
                zeta1 = eval(zeta);
                set(b, 'Value', zeta1);
            end
        else
            zeta1 = eval(get(b2, 'String'));
            set(b, 'Value', zeta1);
        end
        racines = polesSystFrac(mu, omega0, omega1, zeta0, zeta1, alpha);
        racines = racines(imag(racines) >= 0);
        set(points, 'XData', real(racines), 'YData', imag(racines));
    end

    function updateAll()
        mu = eval(get(bmu, 'String'));
        omega1 = omega0*eval(get(bw, 'String'));
        zeta0 = eval(get(bzeta0, 'String'));
        Zeta10 = eval(get(bZeta10, 'String'));
        Zeta11 = eval(get(bZeta11, 'String'));
        Zeta1 = linspace(Zeta10, Zeta11, nZeta1);
        alpha = eval(get(balpha, 'String'));
        
        set(b, 'min', Zeta10, 'max', Zeta11);
        
        racines = nan(4, nZeta1);
        for k = 1:nZeta1
            zeta = Zeta1(k);
            racines(:,k) = polesSystFrac(mu, omega0, omega1, zeta0, zeta, alpha);
        end
        racines = racines(~isnan(racines));
        racines = racines(imag(racines) >= 0); %on garde seulement les poles à partie imaginaire positive pour l'affichage
        set(poles, 'XData', real(racines)', 'YData', imag(racines)');
        
        updateZeta();
    end



end
