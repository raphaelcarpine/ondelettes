function audioPlayer(T, X, initFcn, updateFcn, closeFcn)
%AUDIOPLAYER Summary of this function goes here
%   Detailed explanation goes here

volume = 1;
timeAcc = 1;

if nargin == 0
    Fs = 1000;
    T = (0:10*Fs) / Fs;
    X = randn(2, length(T));
    volume  = 0.1;
end

if nargin <= 2
    initFcn = @(~) 0;
    updateFcn = @(~) 0;
    closeFcn = @() 0;
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
Fs = Fs * timeAcc;

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
Xaudio = Xaudio * volume;

%% initialization

audioP = audioplayer(transpose(Xaudio), Fs);
initFcn(T(1));


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

%
set(tiTxt, 'Units', 'normalized');
set(tfBTxt, 'Units', 'normalized');
set(tTxt, 'Units', 'normalized');
set(timeSlider, 'Units', 'normalized');

%%

kt0 = 0; % index of t0

    function setTime(xt)
        stop(audioP);
        set(playPauseBut, 'String', 'play');
        
        kt0 = floor(xt*length(T));
        if kt0 == length(T)
            kt0 = length(T) - 1;
        end
        
        audioP = audioplayer(transpose(Xaudio(:, kt0+1:end)), Fs);
        initAudioP(audioP);
        
        set(timeSlider, 'Value', xt);
        t = T(1) + xt * (T(end)-T(1));
        set(tTxt, 'String', char(duration(seconds(t), 'Format', 'mm:ss')));
        
        updateFcn(t);
    end

    function updateTime()
        disp(timeSlider.ButtonDownFcn);
        kt = kt0 + get(audioP, 'CurrentSample');
        t = T(kt);
        xt = (t-T(1)) / (T(end)-T(1));
        set(timeSlider, 'Value', xt);
        set(tTxt, 'String', char(duration(seconds(t), 'Format', 'mm:ss')));
        
        updateFcn(t);
    end

    function playPauseFunc()
        flag = isplaying(audioP);
        
        if flag
            set(playPauseBut, 'String', 'play');
            set(audioP, 'UserData', 'paused');
            updateTime();
            pause(audioP);
        else
            set(playPauseBut, 'String', 'pause');
            resume(audioP);
        end
            
    end

    function stopFunc()
        if ~isempty(audioP.UserData)
            set(audioP, 'UserData', []);
            return
        end
        pause(audioP);
        
        kt0 = 0;
        audioP = audioplayer(transpose(Xaudio), Fs);
        initAudioP(audioP);
        
        setTime(0);
        set(playPauseBut, 'String', 'play');
    end

%%

timeSlider.Callback = @(hObject,~) setTime(get(hObject, 'Value'));

playPauseBut.Callback = @(~,~) playPauseFunc();

stopBut.Callback = @(~,~) stopFunc();


    function keyboardShortcuts(~, event)
        switch event.Key
            case 'space'
                playPauseFunc();
        end
    end
set(fig, 'KeyPressFcn', @keyboardShortcuts);

%%

    function initAudioP(audiop)
        set(audiop, 'TimerFcn', @(~,~) updateTime());
%         set(audiop, 'StopFcn', @(~,~) updateTime());
        set(audiop, 'StopFcn', @(~,~) stopFunc());
        set(audiop, 'StartFcn', @(~,~) updateTime());
        % set(audiop, 'TimePeriod', 0.02);
    end

initAudioP(audioP);

%%

    function closeFigFcn(~, ~)
        try
            closeFcn();
        catch
        end
        delete(fig)
    end

fig.CloseRequestFcn = @closeFigFcn;


end