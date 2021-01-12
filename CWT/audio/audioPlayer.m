function audioPlayer(dataPlt, linkedAxes)
%AUDIOPLAYER Summary of this function goes here
%   Detailed explanation goes here

volume = 1;
timeSpeed = 1;

if nargin == 0 % test
    Fs = 1000;
    T = (0:10*Fs) / Fs;
    X = randn(2, length(T));
    volume  = 0.2;
    figure;
    dataPlt = plot(T, X);
end

if nargin <= 1
    linkedAxes = [dataPlt.Parent];
    linkedAxes = unique(linkedAxes);
end

T = [];
X = [];
Xaudio = [];

    function getTX()
        
        try
            T = get(dataPlt, 'XData');
            X = get(dataPlt, 'YData');
            for iT = 1:length(T)
                T{iT} = T{iT}.';
            end
            T = [T{:}].';
            for iX = 1:length(X)
                X{iX} = X{iX}.';
            end
            X = [X{:}].';
            
            T = T(1, :);
        catch
            % deletes plots
        end
        
        %% array size
        
        % array size
        if size(T, 1) == 1 && size(T, 2) == size(X, 2)
            % ok
        elseif size(T, 1) == size(X, 1) && size(T, 2) == size(X, 2)
            dt = mean(diff(T(1, :)));
            for k_t = 2:size(T, 1)
                if max(abs(T(k_t, :) - T(1, :))) > dt * 1e-3
                    error('different time arrays');
                end
            end
            T = T(1, :);
        else
            error(sprintf('array size error (size(t) = [%d, %d] & size(X) = [%d, %d]', [size(T), size(X)]));
        end
        
        % time step
        if any(abs(diff(T)/mean(diff(T)) - 1) > 1e-3)
            error('non-constant time step');
        end
        
        Fs = 1/mean(diff(T));
        
        %% begining
        
        n_begining = round(Fs);
        T = [(-n_begining:-1)/Fs + T(1), T];
        % X = [zeros(size(X, 1), n_begining), X]; % hard begining
        X = [X(:, 1) * (0:n_begining-1)/n_begining, X]; % soft begining
        
        %% stereo
        
        if size(X, 1) == 1
            Xaudio = [X; X];
        elseif size(X, 1) == 2
            Xaudio = X;
        else
            Xaudio = X(1:2, :);
        end
        
        % normalization
        Xaudio = Xaudio - min(Xaudio, [], 'all');
        Xaudio = Xaudio / max(Xaudio, [], 'all');
        Xaudio = 2*Xaudio - 1;
        
    end

%% initialization

getTX()

audioP = audioplayer(volume*transpose(Xaudio), round(timeSpeed*Fs));

audioTimeOnAxes = AudioTimeOnAxes(T(1), linkedAxes);


%% display

L = 100;

fig = figure('Name', 'Audio player', 'numbertitle', 'off');
fig.Units = 'characters';
fig.Position(3) = L;
fig.Position(4) = 7;
fig.MenuBar = 'none';

% time
tiTxt = uicontrol('Parent', fig, 'Style', 'text', 'String', char(duration(seconds(T(1)), 'Format', 'mm:ss')),...
    'Units', 'characters', 'Position', [0, 4.2, 8, 2]);
tfBTxt = uicontrol('Parent', fig, 'Style', 'text', 'String', char(duration(seconds(T(end)), 'Format', 'mm:ss')),...
    'Units', 'characters', 'Position', [L-8, 4.2, 8, 2]);
tTxt = uicontrol('Parent', fig, 'Style', 'text', 'String', char(duration(seconds(T(1)), 'Format', 'mm:ss')),...
    'Units', 'characters', 'Position', [L/2-8/2, 3, 8, 2]);
timeSlider = uicontrol('Parent', fig, 'Style', 'slider',...
    'Units', 'characters', 'Position', [8, 5, L-8*2, 1.3]);

% play button
playPauseBut = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'play',...
    'Units', 'characters', 'Position', [3, 1, 8, 3]);

% stop button
stopBut = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'stop',...
    'Units', 'characters', 'Position', [13, 1, 8, 3]);

% volume
uicontrol('Parent', fig, 'Style', 'text', 'String', 'volume',...
    'Units', 'characters', 'Position', [25, 2.5, 20, 1.5]);
volumeSlider = uicontrol('Parent', fig, 'Style', 'slider', 'Value', volume,...
    'Units', 'characters', 'Position', [25, 1, 20, 1.1]);

% speed
uicontrol('Parent', fig, 'Style', 'text', 'String', 'speed',...
    'Units', 'characters', 'Position', [50, 2, 10, 1]);
speedInput = uicontrol('Parent', fig, 'Style', 'edit', 'String', num2str(timeSpeed),...
    'Units', 'characters', 'Position', [60, 1.5, 6, 2]);

% add linked axes
addAxesInput = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'add linked axes',...
    'Units', 'characters', 'Position', [75, 2.5, 20, 1.8]);
addAxesInput.Callback = @(~,~) audioTimeOnAxes.addAxes(findLinkedAxesAudio());

% sound generator
soundGeneratorInput = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String','sound generator',...
    'Units', 'characters', 'Position', [75, 0.5, 20, 1.8]);
soundGeneratorInput.Callback = @(~,~) generateSoundMenu();

%
set(tiTxt, 'Units', 'normalized');
set(tfBTxt, 'Units', 'normalized');
set(tTxt, 'Units', 'normalized');
set(timeSlider, 'Units', 'normalized');

%%

kt0 = 0; % index of t0


    function initAudioP(wasPlaying)
        if nargin == 0
            wasPlaying = isplaying(audioP);
        end
        
        audioP = audioplayer(volume * transpose(Xaudio(:, kt0+1:end)), round(timeSpeed*Fs));
        
        set(audioP, 'StartFcn', @(~,~) updateTime(true));
        set(audioP, 'TimerFcn', @(~,~) updateTime());
        set(audioP, 'StopFcn', @(~,~) updateTime());
        
        if wasPlaying
            resume(audioP);
        end
    end

    function setTime(xt, wasPlaying)
        if nargin == 1
            wasPlaying = isplaying(audioP);
        end
        
        stop(audioP);
%         set(playPauseBut, 'String', 'play');
        
        kt0 = floor(xt*length(T));
        if kt0 == length(T)
            kt0 = length(T) - 1;
        elseif kt0 == 0
            kt0 = 1;
        end
        
        initAudioP(wasPlaying);
        
        set(timeSlider, 'Value', xt);
        t = T(1) + xt * (T(end)-T(1));
        set(tTxt, 'String', char(duration(seconds(t), 'Format', 'mm:ss')));
        
        audioTimeOnAxes.updateT(t);
    end

    function updateTime(startingFlag)
        if nargin == 0
            startingFlag = false;
        end
        
        kt = kt0 + get(audioP, 'CurrentSample');
        kt = min(kt, length(T));
        t = T(kt);
        xt = (t-T(1)) / (T(end)-T(1));
        set(timeSlider, 'Value', xt);
        set(tTxt, 'String', char(duration(seconds(t), 'Format', 'mm:ss')));
        
        if ~isplaying(audioP) && ~startingFlag
            set(playPauseBut, 'String', 'play');
        end
        
        audioTimeOnAxes.updateT(t);
    end

    function playPauseFunc()
        flag = isplaying(audioP);
        
        if flag
            set(playPauseBut, 'String', 'play');
%             set(audioP, 'UserData', 'paused');
            updateTime();
            pause(audioP);
        else
            set(playPauseBut, 'String', 'pause');
            resume(audioP);
        end
            
    end

    function stopFunc()
%         if ~isempty(get(audioP, 'UserData'))
%             set(audioP, 'UserData', []);
%             return
%         end
        pause(audioP);
        
        getTX()
        
        kt0 = 0;
        
        initAudioP();
        
        setTime(0, false);
        set(playPauseBut, 'String', 'play');
    end

%% callbacks

% time
timeSlider.Callback = @(hObject,~) setTime(get(hObject, 'Value'));

% play/pause
playPauseBut.Callback = @(~,~) playPauseFunc();

% stop
stopBut.Callback = @(~,~) stopFunc();

% volume
    function volumeCallback()
        volume = get(volumeSlider, 'Value');
        kt0 = kt0 + get(audioP, 'CurrentSample') - 1;
%         set(audioP, 'UserData', 'paused');
        initAudioP();
    end
volumeSlider.Callback = @(~,~) volumeCallback();

% speed
    function speedCallback()
        newTimeSpeed = str2double(get(speedInput, 'String'));
        if isnan(newTimeSpeed)
            set(speedInput, 'String', timeSpeed);
        end
        timeSpeed = newTimeSpeed;
        
        kt0 = kt0 + get(audioP, 'CurrentSample') - 1;
%         set(audioP, 'UserData', 'paused');
        initAudioP();
    end
speedInput.Callback = @(~,~) speedCallback();

% shortcuts
    function keyboardShortcuts(~, event)
        switch event.Key
            case 'space'
                playPauseFunc();
        end
    end
set(fig, 'KeyPressFcn', @keyboardShortcuts);
set(tiTxt, 'KeyPressFcn', @keyboardShortcuts);
set(tfBTxt, 'KeyPressFcn', @keyboardShortcuts);
set(tTxt, 'KeyPressFcn', @keyboardShortcuts);
set(timeSlider, 'KeyPressFcn', @keyboardShortcuts);
% set(playPauseBut, 'KeyPressFcn', @keyboardShortcuts);
% set(stopBut, 'KeyPressFcn', @keyboardShortcuts);
set(volumeSlider, 'KeyPressFcn', @keyboardShortcuts);








%% initialization


initAudioP();

%% close acallback

    function closeFigFcn(~, ~)
        try
            audioTimeOnAxes.close();
        catch
        end
        try
            delete(audioP);
        catch
        end
        delete(fig)
    end

fig.CloseRequestFcn = @closeFigFcn;


end