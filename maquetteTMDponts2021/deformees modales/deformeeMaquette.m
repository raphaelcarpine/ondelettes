function [fctDefModale, fctDefModaleAnimation] = deformeeMaquette(positionCapteurs, orientationCapteurs)
%DEFMODALES Summary of this function goes here
%   positionCapteurs = [pos1x, pos2x,... ; pos1z, pos2z,... ]
%   orientationCapteurs = {'xy', 'xyz',... }

saveAsGif = true;

if length(orientationCapteurs) ~= length(positionCapteurs)
    error('length(orientationCapteurs) ~= length(positionCapteurs)');
end

if size(positionCapteurs, 1) == 2
    % position xz => position xyz
    positionCapteurs = [positionCapteurs(1, :); zeros(1, size(positionCapteurs, 2)); positionCapteurs(2, :)];
end

%% dimensions

% maquette
H = 212; % hauteur totale
h = 88; % hauteur de la poutre
l = 302.5; % demie longueur de la poutre


% sol
Dsol = [100 60];

%% capteurs poutre vs pile

Ipile = find((positionCapteurs(1, :) == 0) & (positionCapteurs(3, :) ~= h));
Ipoutre = find((positionCapteurs(3, :) == h) & (positionCapteurs(1, :) ~= 0));
Icentre = find((positionCapteurs(1, :) == 0) & (positionCapteurs(3, :) == h));
if length(Icentre) > 1
    error('plusieurs capteurs au centre');
end

%% fonction deformee

coeffDef = 70;

    function fig = fctDefModale0(animated, deformee, nameFig)
        if nargin < 3
            nameFig = 'deformee modale';
        end
        
        if length(deformee) ~= length(cell2mat(orientationCapteurs))
            error('length(deformee) ~= length(orientationCapteurs)');
        end
        if ~animated && any(imag(deformee) ~= 0)
            error('complex mode shape');
        end
        
        if deformee(1) < 0
            deformee = -deformee;
        end
        
        MovingParts = gobjects(1, 0);
        
        %%
        
        fig = figure('Name', nameFig, 'Color', [1 1 1]);
        fig.UserData = {deformee, positionCapteurs, orientationCapteurs}; % save
        %         fig.Position(3:4) = [250 350];
        ax = axes(fig);
        hold(ax, 'on');
        viewAngle = -35 + 180;
        
        % sol
        cs = 0.5*[1 1 1];
        surf(ax, Dsol(1)*0.5*[-1 -1; 1 1], Dsol(2)*0.5*[-1 1; -1 1], zeros(2), 'FaceColor', cs, 'EdgeCOlor', 'none');
        
        % maquette immobile
        c0 = 0.6*[1 1 1];
        lw0 = 3;
        plot3(ax, [0 0], [0 0], [0 H], 'Color', c0, 'LineWidth', lw0);
        plot3(ax, l*[-1 1], [0 0], h*[1 1], 'Color', c0, 'LineWidth', lw0);
        
        % deformee
        [P2pile, P2poutre, P2capteurs] = getPositions(real(deformee));
        
        % construction deformee
        c1 = 0*[1 1 1];
        lw1 = 3;
        % pile
        for kc = 1:size(P2pile, 2)-1
            MovingParts(end+1) = plot3(ax, [P2pile(1, kc), P2pile(1, kc+1)], [P2pile(2, kc), P2pile(2, kc+1)],...
                [P2pile(3, kc), P2pile(3, kc+1)], 'Color', c1, 'LineWidth', lw1);
        end
        % poutre
        for kc = 1:size(P2poutre, 2)-1
            MovingParts(end+1) = plot3(ax, [P2poutre(1, kc), P2poutre(1, kc+1)], [P2poutre(2, kc), P2poutre(2, kc+1)],...
                [P2poutre(3, kc), P2poutre(3, kc+1)], 'Color', c1, 'LineWidth', lw1);
        end
        % capteurs
        for kc = 1:size(P2capteurs, 2)
            MovingParts(end+1) = plot3(ax, P2capteurs(1, kc), P2capteurs(2, kc), P2capteurs(3, kc), 'o', 'Color', 'red');
        end
        
        % display
        view(ax, viewAngle, 20);
        
        daspect(ax, [1 1 1]);
        set(ax, 'PositionConstraint', 'innerposition');
        set(ax, 'XLim', l*[-1.1 1.1]);
        set(ax, 'YLim', Dsol(2)*[-1 1]);
        set(ax, 'ZLim', H*[0 1.1]);
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
        
        if saveAsGif
            answer = inputdlg('Enter gif saving path', 'GIF saving');
            if ~isempty(answer)
                gifSavingPath = answer{1};
            else
                saveAsGif = false;
            end
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
                    
                    [P2pile, P2poutre, P2capteurs] = getPositions(realShape);
                    
                    % construction deformee
                    km = 1;
                    % pile
                    for kc = 1:size(P2pile, 2)-1
                        MovingParts(km).XData = [P2pile(1, kc), P2pile(1, kc+1)];
                        MovingParts(km).YData = [P2pile(2, kc), P2pile(2, kc+1)];
                        MovingParts(km).ZData = [P2pile(3, kc), P2pile(3, kc+1)];
                        km = km+1;
                    end
                    % poutre
                    for kc = 1:size(P2poutre, 2)-1
                        MovingParts(km).XData = [P2poutre(1, kc), P2poutre(1, kc+1)];
                        MovingParts(km).YData = [P2poutre(2, kc), P2poutre(2, kc+1)];
                        MovingParts(km).ZData = [P2poutre(3, kc), P2poutre(3, kc+1)];
                        km = km+1;
                    end
                    % capteurs
                    for kc = 1:size(P2capteurs, 2)
                        MovingParts(km).XData = P2capteurs(1, kc);
                        MovingParts(km).YData = P2capteurs(2, kc);
                        MovingParts(km).ZData = P2capteurs(3, kc);
                        km = km+1;
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
                
                [P2pile, P2poutre, P2capteurs] = getPositions(realShape);
                
                % construction deformee
                km = 1;
                % pile
                for kc = 1:size(P2pile, 2)-1
                    MovingParts(km).XData = [P2pile(1, kc), P2pile(1, kc+1)];
                    MovingParts(km).YData = [P2pile(2, kc), P2pile(2, kc+1)];
                    MovingParts(km).ZData = [P2pile(3, kc), P2pile(3, kc+1)];
                    km = km+1;
                end
                % poutre
                for kc = 1:size(P2poutre, 2)-1
                    MovingParts(km).XData = [P2poutre(1, kc), P2poutre(1, kc+1)];
                    MovingParts(km).YData = [P2poutre(2, kc), P2poutre(2, kc+1)];
                    MovingParts(km).ZData = [P2poutre(3, kc), P2poutre(3, kc+1)];
                    km = km+1;
                end
                % capteurs
                for kc = 1:size(P2capteurs, 2)
                    MovingParts(km).XData = P2capteurs(1, kc);
                    MovingParts(km).YData = P2capteurs(2, kc);
                    MovingParts(km).ZData = P2capteurs(3, kc);
                    km = km+1;
                end
                drawnow
            catch
            end
        end
        
        if saveAsGif
            gifFunc();
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
        
        function gifFunc
            Nframes = ceil(fps/freq);
            tau = 1/fps;
            t = (0:Nframes-1)/Nframes;
            
            gif(gifSavingPath, 'DelayTime', tau, 'frame', ax, 'nodither');
            
            for k = 1:Nframes
                realShape = real(deformee * exp(2i*pi*t(k)));
                
                [P2pile, P2poutre, P2capteurs] = getPositions(realShape);
                
                % construction deformee
                km = 1;
                % pile
                for kc = 1:size(P2pile, 2)-1
                    MovingParts(km).XData = [P2pile(1, kc), P2pile(1, kc+1)];
                    MovingParts(km).YData = [P2pile(2, kc), P2pile(2, kc+1)];
                    MovingParts(km).ZData = [P2pile(3, kc), P2pile(3, kc+1)];
                    km = km+1;
                end
                % poutre
                for kc = 1:size(P2poutre, 2)-1
                    MovingParts(km).XData = [P2poutre(1, kc), P2poutre(1, kc+1)];
                    MovingParts(km).YData = [P2poutre(2, kc), P2poutre(2, kc+1)];
                    MovingParts(km).ZData = [P2poutre(3, kc), P2poutre(3, kc+1)];
                    km = km+1;
                end
                % capteurs
                for kc = 1:size(P2capteurs, 2)
                    MovingParts(km).XData = P2capteurs(1, kc);
                    MovingParts(km).YData = P2capteurs(2, kc);
                    MovingParts(km).ZData = P2capteurs(3, kc);
                    km = km+1;
                end
                drawnow
                
                gif;
            end
        end
        
    end


    function [P2pile, P2poutre, P2capteurs] = getPositions(deformee0)
        % centre
        Pcentre = [0; 0; h]; % position
        if ~isempty(Icentre)
            Dcentre = zeros(3, 1); % deplacement
            for korient = 1:length(orientationCapteurs{Icentre})
                switch orientationCapteurs{Icentre}(korient)
                    case 'x'
                        Dcentre(1) = deformee0(length(cell2mat(orientationCapteurs(1:Icentre-1))) + korient);
                    case 'y'
                        Dcentre(2) = deformee0(length(cell2mat(orientationCapteurs(1:Icentre-1))) + korient);
                    case 'z'
                        Dcentre(3) = deformee0(length(cell2mat(orientationCapteurs(1:Icentre-1))) + korient);
                    otherwise
                        error('direction non reconnue');
                end
            end
        else
            % recherche capteurs poutres les plus proches
            captProches1 = find((positionCapteurs(1, :) < 0) & (positionCapteurs(3, :) == h));
            [~, I] = max(positionCapteurs(1, captProches1));
            captProches(1) = captProches1(I);
            captProches2 = find((positionCapteurs(1, :) > 0) & (positionCapteurs(3, :) == h));
            [~, I] = min(positionCapteurs(1, captProches2));
            captProches(2) = captProches2(I);
            DcaptProches = zeros(3, 2);
            for kc = 1:2
                for korient = 1:length(orientationCapteurs{captProches(kc)})
                    switch orientationCapteurs{captProches(kc)}(korient)
                        case 'x'
                            DcaptProches(1, kc) = deformee0(length(cell2mat(orientationCapteurs(1:captProches(kc)-1))) + korient);
                        case 'y'
                            DcaptProches(2, kc) = deformee0(length(cell2mat(orientationCapteurs(1:captProches(kc)-1))) + korient);
                            % pas de déplacment selon z au centre a priori
                    end
                end
            end
            CoeffsCaptProches = flip(abs(positionCapteurs(1, captProches)));
            CoeffsCaptProches = CoeffsCaptProches / sum(CoeffsCaptProches);
            
            Dcentre = sum(DcaptProches * CoeffsCaptProches.', 2);
        end
        Dcentre = coeffDef * Dcentre;
        P2centre = Pcentre + Dcentre;
        
        % pile
        Ppile0 = positionCapteurs(:, Ipile);
        Dpile0 = zeros(3, length(Ipile));
        for kc = 1:length(Ipile)
            for korient = 1:length(orientationCapteurs{Ipile(kc)})
                switch orientationCapteurs{Ipile(kc)}(korient)
                    case 'x'
                        Dpile0(1, kc) = deformee0(length(cell2mat(orientationCapteurs(1:Ipile(kc)-1))) + korient);
                    case 'y'
                        Dpile0(2, kc) = deformee0(length(cell2mat(orientationCapteurs(1:Ipile(kc)-1))) + korient);
                    case 'z'
                        Dpile0(3, kc) = deformee0(length(cell2mat(orientationCapteurs(1:Ipile(kc)-1))) + korient);
                    otherwise
                        error('direction non reconnue');
                end
            end
        end
        Dpile0 = coeffDef * Dpile0;
        P2pile0 = Ppile0 + Dpile0;
        
        Ppile = [Ppile0, Pcentre, [0;0;0]];
        Dpile = [Dpile0, Dcentre, [0;0;0]];
        [~, I] = sort(Ppile(3, :));
        Ppile = Ppile(:, I);
        Dpile = Dpile(:, I);
        P2pile = Ppile + Dpile;
        
        % poutre
        Ppoutre0 = positionCapteurs(:, Ipoutre);
        Dpoutre0 = zeros(3, length(Ipoutre));
        for kc = 1:length(Ipoutre)
            for korient = 1:length(orientationCapteurs{Ipoutre(kc)})
                switch orientationCapteurs{Ipoutre(kc)}(korient)
                    case 'x'
                        Dpoutre0(1, kc) = deformee0(length(cell2mat(orientationCapteurs(1:Ipoutre(kc)-1))) + korient);
                    case 'y'
                        Dpoutre0(2, kc) = deformee0(length(cell2mat(orientationCapteurs(1:Ipoutre(kc)-1))) + korient);
                    case 'z'
                        Dpoutre0(3, kc) = deformee0(length(cell2mat(orientationCapteurs(1:Ipoutre(kc)-1))) + korient);
                    otherwise
                        error('direction non reconnue');
                end
            end
        end
        Dpoutre0 = coeffDef * Dpoutre0;
        P2poutre0 = Ppoutre0 + Dpoutre0;
        
        Ppoutre = [Ppoutre0, Pcentre];
        Dpoutre = [Dpoutre0, Dcentre];
        [~, I] = sort(Ppoutre(1, :));
        Ppoutre = Ppoutre(:, I);
        Dpoutre = Dpoutre(:, I);
        P2poutre = Ppoutre + Dpoutre;
        
        % capteurs
        P2capteurs = [P2pile0, P2poutre0];
        if ~isempty(Icentre)
            P2capteurs = [P2capteurs, P2centre];
        end
    end


fctDefModale = @(varargin) fctDefModale0(false, varargin{:});
fctDefModaleAnimation = @(varargin) fctDefModale0(true, varargin{:});


end

