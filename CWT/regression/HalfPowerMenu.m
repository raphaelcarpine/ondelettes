function [fmax, zeta] = HalfPowerMenu()
%WaveletMenu Summary of this function goes here
%   Detailed explanation goes here

fmax = nan;
zeta = nan;

quadraticMax = true;

quadraticFFT = false;

%% fig

% lignes pans
Hpans = [4, 4, 10];

fig = figure('Name', 'Half power menu', 'numbertitle', 'off', 'Resize', 'off');
fig.Units = 'characters';
fig.Position(3) = 50;
fig.Position(4) = sum(Hpans);
fig.Position(2) = fig.Position(2) - fig.Position(4) - 2;
fig.MenuBar = 'none';



%% bouton regression et panneaux param et sorties

linePan = uipanel('Parent', fig, 'Units', 'normalized');
maxPan = uipanel('Parent', fig, 'Units', 'normalized');
resultsPan = uipanel('Parent', fig, 'Units', 'normalized');

margin = 0.02;
hpans = (1 - margin*(length(Hpans)+1)) * Hpans/sum(Hpans);
zpans = margin * ones(size(hpans));
for iz = length(zpans)-1:-1:1
    zpans(iz) = zpans(iz+1) + hpans(iz+1) + margin;
end

linePan.Position = [0.02, zpans(1), 0.96, hpans(1)];
maxPan.Position = [0.02, zpans(2), 0.96, hpans(2)];
resultsPan.Position = [0.02, zpans(3), 0.96, hpans(3)];



%% ligne � fitter

line = [];

lines = [];
kline = 0;

highlighted = [];
lineWidth = [];


lineSelect = uicontrol('Parent', linePan, 'Units', 'normalized','Style','togglebutton', 'String', 'select line');
prevBut = uicontrol('Parent', linePan, 'Units', 'normalized','Style','pushbutton', 'String', '<');
nextBut = uicontrol('Parent', linePan, 'Units', 'normalized','Style','pushbutton', 'String', '>');

lineSelect.Position = [0.01, 0.02, 0.48, 0.96];
prevBut.Position = [0.6, 0.1, 0.15, 0.8];
nextBut.Position = [0.75, 0.1, 0.15, 0.8];


    function highlightLine()
        try
            set(highlighted, 'LineWidth', lineWidth);
            lineWidth = get(line, 'LineWidth');
            highlighted = line;
            set(highlighted, 'LineWidth', 3*lineWidth);
        catch
        end
    end

    function selectFunction(selecting)
        if selecting
            lines = findobj('Type', 'line');
            kline = 1;
            line = lines(1:min(1,end));
        else
            line = [];
            lines = [];
        end
        highlightLine();
        
        maxSelect.Value = false;
        selectFunctionMax(maxSelect.Value);
    end

    function nextprev(next)
        if ~isempty(lines)
            kline = mod(kline + next -1, length(lines)) + 1;
            line = lines(kline);
            highlightLine();
        end
    end

lineSelect.Callback = @(~,~) selectFunction(lineSelect.Value);
prevBut.Callback = @(~,~) nextprev(-1);
nextBut.Callback = @(~,~) nextprev(1);

    function closeReg()
        try
            selectFunction(false);
        catch
        end
        delete(fig);
    end

fig.CloseRequestFcn = @(~,~) closeReg();



%% max � fitter

localMax = [];

localMaxs = [];
kmax = 0;

ax = 0;

maxMarker = gobjects(1);
halfMarker = gobjects(1);


maxSelect = uicontrol('Parent', maxPan, 'Units', 'normalized','Style','togglebutton', 'String', 'select local maximum');
prevButmax = uicontrol('Parent', maxPan, 'Units', 'normalized','Style','pushbutton', 'String', '<');
nextButmax = uicontrol('Parent', maxPan, 'Units', 'normalized','Style','pushbutton', 'String', '>');

maxSelect.Position = [0.01, 0.02, 0.48, 0.96];
prevButmax.Position = [0.6, 0.1, 0.15, 0.8];
nextButmax.Position = [0.75, 0.1, 0.15, 0.8];


    function highlightMax()
        delete(maxMarker);
        delete(halfMarker);
        fmax = nan;
        zeta = nan;
        displayFmaxZeta();
        if ~isempty(localMax)
            computeHalfPower();
        end
    end

    function selectFunctionMax(selecting)
        if selecting && updateXYAxes()
            localMaxs = findLocalMaxs(Y);
            kmax = 1;
            localMax = localMaxs(1:min(1,end)); % element #1 if length ~= 0, [] otherwise
        else
            localMax = [];
            localMaxs = [];
        end
        highlightMax();
    end

    function nextprevMax(next)
        if ~isempty(lines)
            kmax = mod(kmax + next -1, length(localMaxs)) + 1;
            localMax = localMaxs(kmax);
            highlightMax();
        end
    end

    function Lm = findLocalMaxs(L)
        Lm = [];
        for k = 2:length(L)-1
            if L(k-1) < L(k) && L(k) >= L(k+1)
                Lm(end+1) = k;
            end
        end
        [~, I] = sort(L(Lm), 'descend');
        Lm = Lm(I);
    end

    function computeHalfPower()
        % fmax
        if quadraticMax
            [fmax, Ymax] = localMax3Points(X(localMaxs(kmax)-1:localMaxs(kmax)+1).', Y(localMaxs(kmax)-1:localMaxs(kmax)+1).');
        else
            fmax = X(localMaxs(kmax));
            Ymax = Y(localMaxs(kmax));
        end
        
        % zeta
        if quadraticFFT
            H = Ymax/2;
        else
            H = Ymax/sqrt(2);
        end
        k1 = localMaxs(kmax)-1;
        while k1 > 0 && Y(k1) > H
            k1 = k1 - 1;
        end
        k2 = localMaxs(kmax)+1;
        while k2 <= length(Y) && Y(k2) > H
            k2 = k2 + 1;
        end
        if k1 < 1 || k2 > length(Y)
            zeta = nan;
        else
            fh1 = X(k1) + (H-Y(k1))*(X(k1+1)-X(k1))/(Y(k1+1)-Y(k1));
            fh2 = X(k2-1) + (H-Y(k2-1))*(X(k2)-X(k2-1))/(Y(k2)-Y(k2-1));
            Q = fmax/(fh2-fh1);
            zeta = 1/(2*Q);
        end
        
        % display
        displayFmaxZeta()
        
        % plot
        maxMarker = scatter(ax, fmax, Ymax, '+', 'MarkerEdgeColor', 'red', 'LineWidth', 2);
        halfMarker = yline(ax, H, 'Color', 0.5*[1 1 1]);
    end

maxSelect.Callback = @(~,~) selectFunctionMax(maxSelect.Value);
prevButmax.Callback = @(~,~) nextprevMax(-1);
nextButmax.Callback = @(~,~) nextprevMax(1);


%% quadratic fft

bg = uibuttongroup('Parent', resultsPan, 'Units', 'characters',...
    'Position', [1 6.5 22 2.3], 'SelectionChangedFcn', @bselection);

butFFT = uicontrol(bg, 'Units', 'characters',...
    'Style', 'radiobutton', 'String', 'FFT');
butFFT2 = uicontrol(bg, 'Units', 'characters',...
    'Style', 'radiobutton', 'String', 'FFT�', 'Value', quadraticFFT);

butFFT.Position = [1 0.5 10 1];
butFFT2.Position = [12 0.5 10 1];

    function bselection(~, ~)
        if butFFT2.Value
            quadraticFFT = true;
        else
            quadraticFFT = false;
        end
        highlightMax();
    end


%% display

FmaxDisp = uicontrol(resultsPan, 'Style', 'text', 'Units', 'characters',...
    'FontSize', 10, 'HorizontalAlignment', 'left');
ZetaDisp = uicontrol(resultsPan, 'Style', 'text', 'Units', 'characters',...
    'FontSize', 10, 'HorizontalAlignment', 'left');

FmaxDisp.Position = [2, 4, 40, 1.5];
ZetaDisp.Position = [2, 2, 40, 1.5];

    function displayFmaxZeta()
        FmaxDisp.String = sprintf('Frequency: %.3f Hz', fmax);
        ZetaDisp.String = sprintf('Damping: %.2f %%', 100*zeta);
    end

displayFmaxZeta();

%% donnees

X = nan;
Y = nan;

    function ok = updateXYAxes()
        %         line = findobj('Type', 'line');
        if ~isempty(line)
            %             line = line(1);
            X = get(line, 'XData');
            Y = get(line, 'YData');
            
            X = X(~isnan(Y)); % on enl�ve les valeurs inappropri�es
            Y = Y(~isnan(Y));
            X = X(~isnan(X));
            Y = Y(~isnan(X));
            
            ax = get(line, 'Parent');
            hold(ax, 'on');
            ok = ~isempty(X);
        else
            ok = false;
        end
    end





%% raccourcis

%     function keyboardShortcuts(~, event)
%         if isequal(event.Key, 'return') || isequal(event.Key, 'space')
%             computeReg('reg');
%         elseif isequal(event.Key, 'shift') || isequal(event.Key, 's')
%             lineSelect.Value = ~ lineSelect.Value;
%             selectFunction(lineSelect.Value);
%         elseif isequal(event.Key, 'p')
%             plotFunc();
%         elseif isequal(event.Key, 'd')
%             deletePlots();
%         elseif isequal(event.Key, 'leftarrow')
%             nextprev(-1);
%         elseif isequal(event.Key, 'rightarrow')
%             nextprev(1);
%         elseif isequal(event.Key, 'delete') || isequal(event.Key, 'backspace')
%             close(fig);
%         end
%     end
% set(fig, 'KeyPressFcn', @keyboardShortcuts);


%% fermeture

if nargout > 0
    waitfor(fig);
end


end