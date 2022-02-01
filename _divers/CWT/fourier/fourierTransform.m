function [f, Fx] = fourierTransform(t, x, varargin)
%FOURIERTRANSFORM Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
%% default param
AveragingDef = false;
AveragingNbDef = 10;
WindowDef = 'rectangular'; % 'rectangular', 'hamming', 'exponential'
WindowParamsDef = [];
WindowSequenceDef = 'symetric'; % 'symetric', 'periodic'

addRequired(p,'t', @isrow);
addRequired(p,'x', @isrow);
addParameter(p, 'Averaging', AveragingDef);
addParameter(p, 'AveragingNb', AveragingNbDef);
addParameter(p, 'Window', WindowDef);
addParameter(p, 'WindowParams', WindowParamsDef);
addParameter(p, 'WindowSequence', WindowSequenceDef);

parse(p, t, x, varargin{:});

Averaging = p.Results.Averaging;
AveragingNb = p.Results.AveragingNb;
Window = p.Results.Window;
WindowParams = p.Results.WindowParams;
WindowSequence = p.Results.WindowSequence;

%%
N = size(t, 2);
dt = mean(diff(t));
if any(abs(diff(t)/dt - 1) > 1e-3)
    error('non-constant time step');
end


%% fft

if ~Averaging
    % window
    x = windowX(x);
    
    % fft
    Fx = fft(x) / N;
    Fx = Fx(1:floor(end/2));
    Fx(2:end) = 2*Fx(2:end);
    f = 1/(N*dt) * (0:length(Fx)-1);
else
    % cut
    n = floor(N/AveragingNb);
    x = x(1:end-mod(N, AveragingNb));
    x = transpose(reshape(x, [n, AveragingNb]));
    
    % window
    x = windowX(x);
    
    Fx = transpose(fft(transpose(x)));
    Fx = abs(Fx).^2;
    Fx = sqrt(mean(Fx, 1));
    
    Fx = Fx / n;
    Fx = Fx(1:floor(end/2));
    Fx(2:end) = 2*Fx(2:end);
    f = 1/(n*dt) * (0:length(Fx)-1);
end

% windowing

function x = windowX(x)
    Nx = size(x, 2);
    
    switch WindowSequence
        case 'symetric'
            Nw = Nx+1;
        case 'periodic'
            Nw = Nx;
    end
    
    switch Window
        case 'rectangular'
            w = ones(1, Nx);
        case 'hamming'
            a0 = 25/46;
            w = a0 - (1-a0)*cos(2*pi*(0:Nx-1)/Nw);
            w = w / mean(w);
        case 'exponential'
            tau = WindowParams / dt;
            w = exp(-(0:Nx-1)/tau);
    end
    
    x = x.*w;
end


end
