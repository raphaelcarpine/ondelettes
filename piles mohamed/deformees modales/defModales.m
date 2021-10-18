function [fctDefModale, fctDefModaleAnimation] = defModales(pile, chNames, nameFig, zSable)
%DEFMODALES Summary of this function goes here
%   pile : 1-10
%   chNames = {'acc1x', 'acc1y', 'acc1z', 'acc2x', 'acc2y', 'acc3x',... }

if nargin < 3
    nameFig = sprintf('deformée pile %u', pile);
end
if nargin < 4
    zSable = nan;
end

%% dimensions

% pile
Dx = [16, 16, 10, 10, 10, 16, 16, 10, 10, 10];
Dy = [16, 16, 15, 15, 30, 16, 16, 15, 15, 30];
Dz = [75, 100, 75, 100, 100, 100, 75, 100, 75, 100];
cyl = [true, true, false, false, false, true, true, false, false, false]; % cylindre

% semelle
Dx_sem = nan(size(Dx));
Dy_sem = nan(size(Dy));
Dz_sem = zeros(size(Dz));
for kp = 6:length(Dx)
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
Dx = Dx(pile);
Dy = Dy(pile);
Dz = Dz(pile);
D = [Dx, Dy, Dz];
cyl = cyl(pile);
% semelle
Dx_sem = Dx_sem(pile);
Dy_sem = Dy_sem(pile);
Dz_sem = Dz_sem(pile);
D_sem = [Dx_sem, Dy_sem, Dz_sem];

% capteurs
if Dz == 75
    Cx = [0, 1, 1, 1] * Dx/2;
    Cy = [0 0 0 0];
    Cz = Dz + Dz_sem + [0, -5, -20, -35];
elseif Dz == 100
    Cx = [0, 1, 1, 1, 1] * Dx/2;
    Cy = [0 0 0 0 0];
    Cz = Dz + Dz_sem + [0, -10, -30, -50, -70];
end
if cyl
    Cx(1) = -Dx/2;
else
    Cy(1) = -Dy/2;
end
C = [Cx; Cy; Cz];

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

% ch 2x
Kch2x = nan;
for kch = 1:length(chNames)
    if strcmp(chNames{kch}(end-1:end), '1x')
        Kch2x = kch;
        break
    end
end
if isnan(Kch2x)
    for kch = 1:length(chNames)
        if strcmp(chNames{kch}(end-1:end), '1y')
            Kch2x = kch;
            break
        end
    end
end
if isnan(Kch2x)
    error(' ');
end

%% fonction deformee

coeffArrows = 40;

    function fig = fctDefModale0(animated, deformee, nameFig2)
        if nargin < 3
            nameFig2 = nameFig;
        end
        
        if length(deformee) ~= length(chNames)
            error('length(deformee) ~= length(chNames)');
        end
        if ~animated && any(imag(deformee) ~= 0)
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
        viewAngle = -65 + 180;
        
        % pile
        faceColor = 0.7*[1 1 1];
        faceAlpha = 0.65;
        edgeColor = 'black';
        if cyl % cylindre
            th = linspace(0, 2*pi, 51);
            Xcircle = D(1)/2 * cos(th);
            Ycircle = D(1)/2 * sin(th);
            Zcircle = ones(size(Xcircle));
            fill3(ax, Xcircle, Ycircle, D_sem(3) + 0*Zcircle, faceColor, 'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha);
            fill3(ax, Xcircle, Ycircle, D_sem(3) + D(3)*Zcircle, faceColor, 'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha);
            Xcircle = [1; 1]*Xcircle;
            Ycircle = [1; 1]*Ycircle;
            Zcircle = D_sem(3) + D(3)*[0; 1]*Zcircle;
            surf(ax, Xcircle, Ycircle, Zcircle, 'EdgeColor', 'none', 'FaceColor', faceColor, 'FaceAlpha', faceAlpha);
            plot3(ax, D(1)/2*cos(viewAngle*pi/180)*[1 1], D(1)/2*sin(viewAngle*pi/180)*[1 1], D_sem(3) + D(3)*[0; 1], edgeColor);
            plot3(ax, D(1)/2*cos(viewAngle*pi/180+pi)*[1 1], D(1)/2*sin(viewAngle*pi/180+pi)*[1 1], D_sem(3) + D(3)*[0; 1], edgeColor);
        else
            Xplate = [-0.5, -0.5; -0.5, -0.5; 0.5, 0.5; 0.5, 0.5] * D(1);
            Yplate = [-0.5, 0.5; -0.5, 0.5; -0.5, 0.5; -0.5, 0.5] * D(2);
            Zplate = D_sem(3) + [0 0; 1 1; 1 1; 0 0] * D(3);
            surf(ax, Xplate, Yplate, Zplate, 'EdgeColor', edgeColor, 'FaceColor', faceColor, 'FaceAlpha', faceAlpha);
            Xplate = [-0.5, 0.5; -0.5, 0.5; -0.5, 0.5; -0.5, 0.5] * D(1);
            Yplate = [-0.5, -0.5; -0.5, -0.5; 0.5, 0.5; 0.5, 0.5] * D(2);
            Zplate = D_sem(3) + [1 1; 0 0; 0 0; 1 1] * D(3);
            surf(ax, Xplate, Yplate, Zplate, 'EdgeColor', edgeColor, 'FaceColor', faceColor, 'FaceAlpha', faceAlpha);
        end
        
        % semelle
        if ~isnan(D_sem(1))
            Xplate = [-0.5, -0.5; -0.5, -0.5; 0.5, 0.5; 0.5, 0.5] * D_sem(1);
            Yplate = [-0.5, 0.5; -0.5, 0.5; -0.5, 0.5; -0.5, 0.5] * D_sem(2);
            Zplate = [0 0; 1 1; 1 1; 0 0] * D_sem(3);
            surf(ax, Xplate, Yplate, Zplate, 'EdgeColor', edgeColor, 'FaceColor', faceColor, 'FaceAlpha', faceAlpha);
            Xplate = [-0.5, 0.5; -0.5, 0.5; -0.5, 0.5; -0.5, 0.5] * D_sem(1);
            Yplate = [-0.5, -0.5; -0.5, -0.5; 0.5, 0.5; 0.5, 0.5] * D_sem(2);
            Zplate = [1 1; 0 0; 0 0; 1 1] * D_sem(3);
            surf(ax, Xplate, Yplate, Zplate, 'EdgeColor', edgeColor, 'FaceColor', faceColor, 'FaceAlpha', faceAlpha);
        end
        
        % sable
        if ~isnan(zSable)
            Dsable = 80;
            Xplate = [-0.5, 0.5; -0.5, 0.5] * Dsable;
            Yplate = [-0.5, -0.5; 0.5, 0.5] * Dsable;
            Zplate = D_sem(3) + zSable*[1 1; 1 1];
            surf(ax, Xplate, Yplate, Zplate, 'EdgeColor', 'none', 'FaceColor', 'yellow', 'FaceAlpha', 0.3);
        end
        
        % capteurs
        accColor = 'black';
        accDim = 2;
        for kcapt1 = 1:size(C, 2)
            XYZ = C(:, kcapt1);
            Xplate = XYZ(1) + ([0 0; 0 0; 1 1; 1 1] - 0.5) * accDim;
            Yplate = XYZ(2) + ([0 1; 0 1; 0 1; 0 1] - 0.5) * accDim;
            Zplate = XYZ(3) + ([0 0; 1 1; 1 1; 0 0] - 0.5) * accDim;
            surf(ax, Xplate, Yplate, Zplate, 'EdgeColor', accColor, 'FaceColor', accColor, 'FaceAlpha', 1);
            Xplate = XYZ(1) + ([0 1; 0 1; 0 1; 0 1] - 0.5) * accDim;
            Yplate = XYZ(2) + ([0 0; 0 0; 1 1; 1 1] - 0.5) * accDim;
            Zplate = XYZ(3) + ([1 1; 0 0; 0 0; 1 1] - 0.5) * accDim;
            surf(ax, Xplate, Yplate, Zplate, 'EdgeColor', accColor, 'FaceColor', accColor, 'FaceAlpha', 1);
        end
        
        % fleches
        arrowColorXYZ = {'red', [.5 0 .5], 'blue'};
        arrowMarker = {'^', 'v'; '>', '<'; '^', 'v'};
        arrowPlots = gobjects(0);
        for kchan = 1:length(deformee)
            kcapt = str2double(chNames{kchan}(end-1));
            XYZ1 = C(:, kcapt);
            XYZ2 = C(:, kcapt) + real(deformee(kchan))*coeffArrows*Dir_C(:, kchan);
            colArrow = arrowColorXYZ{logical(Dir_C(:, kchan))};
            arrowPlots(end+1) = plot3(ax, [XYZ1(1), XYZ2(1)], [XYZ1(2), XYZ2(2)], [XYZ1(3), XYZ2(3)],...
                '-', 'Color', colArrow, 'LineWidth', 2, 'MarkerSize', 2,...
                'Marker', arrowMarker{logical(Dir_C(:, kchan)), (real(deformee(kchan)) < 0) + 1},...
                'MarkerIndices', 2, 'MarkerEdgeCOlor', colArrow, 'MarkerFaceColor', colArrow);
        end
        
        % display
        view(ax, viewAngle, 20);
        
        daspect(ax, [1 1 1]);
        set(ax, 'PositionConstraint', 'innerposition');
        set(ax, 'XLim', 80*[-0.5, 0.5]);
        set(ax, 'YLim', 80*[-0.5, 0.5]);
        set(ax, 'ZLim', [0, D(3) + D_sem(3) + coeffArrows/2]);
        set(ax, 'Clipping', 'off');
        % pbaspect(ax, [plateDim, 0.3*max(plateDim)]);
        set(ax,'visible','off');
        set(ax, 'InnerPosition', [0.04 0.04 0.92 0.92]);
        set(ax, 'ClippingStyle', 'rectangle');
        
%         lightangle(ax, -100, 40);
%         lighting flat
        
        %% animation
        if ~animated
            return
        end
        
        freq = 0.6; % display frequency
        fps = 25; % nb of frames
        
        function animation
            Nframes = ceil(fps/freq);
            tau = 1/fps;
            t = (0:Nframes-1)/Nframes;
            tic;
            while fig == get(groot, 'CurrentFigure')
                try
                    k = floor(toc/tau) +1;
                    realShape = real(deformee * exp(2i*pi*t(mod(k-1, Nframes) +1)));
                    for kchan2 = 1:length(realShape)
                        kcapt = str2double(chNames{kchan2}(end-1));
                        XYZ1 = C(:, kcapt);
                        XYZ2 = C(:, kcapt) + realShape(kchan2)*coeffArrows*Dir_C(:, kchan2);
                        arrowPlots(kchan2).XData = [XYZ1(1), XYZ2(1)];
                        arrowPlots(kchan2).YData = [XYZ1(2), XYZ2(2)];
                        arrowPlots(kchan2).ZData = [XYZ1(3), XYZ2(3)];
                        arrowPlots(kchan2).Marker = arrowMarker{logical(Dir_C(:, kchan2)), (real(realShape(kchan2)) < 0) + 1};
                    end
                    drawnow
                    pause((k+1)*tau - toc);
                catch
                    break
                end
            end
            
            % end
            try
                realShape = real(deformee);
                for kchan2 = 1:length(realShape)
                    kcapt = str2double(chNames{kchan2}(end-1));
                    XYZ1 = C(:, kcapt);
                    XYZ2 = C(:, kcapt) + realShape(kchan2)*coeffArrows*Dir_C(:, kchan2);
                    arrowPlots(kchan2).XData = [XYZ1(1), XYZ2(1)];
                    arrowPlots(kchan2).YData = [XYZ1(2), XYZ2(2)];
                    arrowPlots(kchan2).ZData = [XYZ1(3), XYZ2(3)];
                    arrowPlots(kchan2).Marker = arrowMarker{logical(Dir_C(:, kchan2)), (real(realShape(kchan2)) < 0) + 1};
                end
                drawnow
            catch
            end
        end
        
        clickTxt = annotation(fig, 'textbox', [0 0 1 1], 'String', 'Click to start animation', 'FitBoxToText', 'off',...
            'FontSize', 16, 'BackgroundColor', 'white', 'FaceAlpha', 0.5, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        set(get(ax,'Children'),'HitTest','off')
        
        
        function clickFcn(~,~)
            delete(clickTxt);
            animation;
        end
        
        fig.WindowButtonDownFcn = @clickFcn;
        fig.Interruptible = 'off';
        
    end


fctDefModale = @(varargin) fctDefModale0(false, varargin{:});
fctDefModaleAnimation = @(varargin) fctDefModale0(true, varargin{:});


end

