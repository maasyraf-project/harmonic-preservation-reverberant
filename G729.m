function y = G729( x, Fs, Nw, Nsh )
%G729 Summary of this function goes here
%   Detailed explanation goes here
% frameSize = Nw-Nsh;
VADPar = InitVADPar(Fs,Nw,Nsh);
numFrame  = floor((length(x)-Nw)/Nsh+1);

y = zeros(size(x));

x_start   = 1;
x_end     = Nsh; % the window shift logic is already supported by G729 implementation

for j=1:numFrame
    x_frame = x(x_start:x_end);
    if length(x_frame) < Nsh
        x_frame = [x_frame zeros(1, Nsh - length(x_frame))];
    end
    [v,VADPar] = VAD(x_frame, VADPar);
    y(x_start:x_end) = v;
    x_start = x_start+Nsh;
    x_end   = min(x_end+Nsh,length(x));
end