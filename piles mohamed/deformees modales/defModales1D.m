function fctDefModale = defModales1D(pile, chNames, direction, nameFig, zSable)
%DEFMODALES Summary of this function goes here
%   pile : 1-10
%   chNames = {'acc1x', 'acc1y', 'acc1z', 'acc2x', 'acc2y', 'acc3x',... }

if nargin < 4
    nameFig = sprintf('deformée pile %u', pile);
end
if nargin < 5
    zSable = nan;
end

regLin = true;
plotSemelleDef = true;
rotateSemelleDef = true;

%% dimensions

% pile
Dz = [75, 100, 75, 100, 100, 100, 75, 100, 75, 100];

% semelle
Dx_sem = nan(size(Dz));
Dy_sem = nan(size(Dz));
Dz_sem = zeros(size(Dz));
for kp = 6:length(Dz)
    Dx_sem(kp) = 45;
    Dy_sem(kp) = 45;
    if Dz(kp) == 75
        Dz_sem(kp) = 20;
    elseif Dz(kp) == 100
        Dz_sem(kp) = 30;
    end
end

%% dimensions pile

% pile
Dz = Dz(pile);
% semelle
Dx_sem = Dx_sem(pile);
Dy_sem = Dy_sem(pile);
Dz_sem = Dz_sem(pile);
D_sem = [Dx_sem, Dy_sem, Dz_sem];

% capteurs
if Dz == 75
    Cz = Dz + Dz_sem + [0, -5, -20, -35];
elseif Dz == 100
    Cz = Dz + Dz_sem + [0, -10, -30, -50, -70];
end

% direction capteurs
Dir_C = zeros(3, length(chNames));
for kch = 1:length(chNames)
    if chNames{kch}(end) == 'x'
        Dir_C(1, kch) = 1;
    elseif chNames{kch}(end) == 'y'
        Dir_C(2, kch) = 1;
    elseif chNames{kch}(end) == 'z'
        Dir_C(3, kch) = 1;
    end
end

% ch 1 pour orientation deformee
Kch2x = nan;
for kch = 1:length(chNames)
    if strcmp(chNames{kch}(end-1:end), ['1', direction])
        Kch2x = kch;
        break
    end
end
if isnan(Kch2x)
    error('');
end

%% fonction deformee

coeffArrows = 40;

    function fig = fctDefModale0(deformee, nameFig2)
        if nargin < 3
            nameFig2 = nameFig;
        end
        
        if length(deformee) ~= length(chNames)
            error('length(deformee) ~= length(chNames)');
        end
        if any(imag(deformee) ~= 0)
            error('complex mode shape');
        end
        
        if deformee(Kch2x) < 0
            deformee = -deformee;
        end
        
        %%
        
        fig = figure('Name', nameFig2);
        fig.UserData = {pile, chNames, nameFig2, zSable, deformee}; % save
        fig.Position(3:4) = [250 350];
        ax = axes(fig);
        hold(ax, 'on');
        
        % pile
        cPile = 0.6*[1 1 1];
        plot(ax, [0 0], Dz_sem + [0, Dz], 'Color', cPile, 'LineWidth', 2);
        
        % semelle
        if ~isnan(D_sem(1))
            Xsem = 0.5*[-1 -1 1 1 -1] * Dx_sem;
            Ysem = [0 1 1 0 0] * Dz_sem;
            plot(ax, Xsem, Ysem, 'Color', cPile, 'LineWidth', 2);
        end
        
        % sable
        if ~isnan(zSable)
            Dsable = 80;
            plotSable = plot(ax, 0.5*[-1 1]*Dsable, zSable + Dz_sem + [0 0], 'Color', 0.8*[1 1 0], 'LineWidth', 1);
            uistack(plotSable, 'bottom');
        end
        
        % repere
        cRep = 0.3*[1 1 1];
        xRep = -40;
        zRep = Dz_sem + Dz + 10;
        lRep = 15;
        plot(ax, xRep + [0 0], zRep + [0 lRep], '-^', 'Color', cRep,...
            'MarkerIndices', 2, 'MarkerSize', 3, 'MarkerFaceColor', cRep);
        text(ax, xRep, zRep + lRep, {'z', '', ''}, 'FontSize', 9, 'Color', cRep,...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        plot(ax, xRep + [0 lRep], zRep + [0 0], '->', 'Color', cRep,...
            'MarkerIndices', 2, 'MarkerSize', 3, 'MarkerFaceColor', cRep);
        text(ax, xRep + lRep, zRep, ['   ', direction], 'FontSize', 9, 'Color', cRep,...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        
        % regression deformee
        if regLin
            Zcapt = nan(length(Cz), 1);
            Xcapt = nan(length(Cz), 1);
            for kchan = 1:length(deformee)
                if chNames{kchan}(end) ~= direction
                    continue
                end
                kcapt = str2double(chNames{kchan}(end-1));
                Zcapt(kcapt) = Cz(kcapt);
                Xcapt(kcapt) = real(deformee(kchan))*coeffArrows;
            end
            regLinCoeffs = [ones(size(Zcapt)), Zcapt] \ Xcapt;
            
            cPile = 0.*[1 1 1];
            plot(ax, regLinCoeffs(1) + regLinCoeffs(2)*(Dz_sem + [0, Dz]), Dz_sem + [0, Dz],...
                'Color', cPile, 'LineWidth', 2);
        end
        
        % deformee
        colSensor = 'red';
        sensorPlots = gobjects(0);
        for kchan = 1:length(deformee)
            if chNames{kchan}(end) ~= direction
                continue
            end
            kcapt = str2double(chNames{kchan}(end-1));
            sensorPlots(end+1) = scatter(ax, real(deformee(kchan))*coeffArrows,...
                Cz(kcapt), 50, colSensor, '+', 'LineWidth', 2);
        end
        
        % semelle deformee
        if plotSemelleDef
            if ~isnan(D_sem(1))
                if ~rotateSemelleDef
                    Xsem = 0.5*[-1 -1 1 1 -1] * Dx_sem;
                    Ysem = [0 1 1 0 0] * Dz_sem;
                    Xsem = Xsem + regLinCoeffs(1) + regLinCoeffs(2)*Ysem;
                    plot(ax, Xsem, Ysem, 'Color', cPile, 'LineWidth', 2);
                else
                    c2 = regLinCoeffs(2);
                    c2 = c2/sqrt(1+c2^2); % sin(arctan)
                    c1 = sqrt(1-c2^2);
                    Mrot = [c1 c2; -c2 c1];
                    Xsem = 0.5*[-1 -1 1 1 -1] * Dx_sem;
                    Ysem = [0 1 1 0 0] * Dz_sem;
                    XYsem = Mrot * [Xsem; Ysem - Dz_sem];
                    Xsem = XYsem(1, :);
                    Ysem = XYsem(2, :);
                    Xsem = Xsem + regLinCoeffs(1) + regLinCoeffs(2) * Dz_sem;
                    Ysem = Ysem + Dz_sem;
                    plot(ax, Xsem, Ysem, 'Color', cPile, 'LineWidth', 2);
                end
            end
        end
        
        % display
        daspect(ax, [1 1 1]);
%         set(ax, 'PositionConstraint', 'innerposition');
        set(ax, 'XLim', Dsable*[-0.5, 0.5]);
        if plotSemelleDef && ~isnan(D_sem(1)) && rotateSemelleDef
            YlimMin = -10;
        else
            YlimMin = 0;
        end
        set(ax, 'YLim', [YlimMin, Dz + Dz_sem + 30]);
%         set(ax, 'Clipping', 'off');
        % pbaspect(ax, [plateDim, 0.3*max(plateDim)]);
        set(ax,'visible','off');
        set(ax, 'Clipping', 'off');
        set(ax, 'InnerPosition', [0.04 0.04 0.92 0.92]);
        
        
    end


fctDefModale = @(varargin) fctDefModale0(varargin{:});


end

