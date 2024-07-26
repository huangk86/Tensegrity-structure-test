function [rawSpectrum] = WindowedFFT (InputSignal,WindowName)
% Direct FFT with rectangular window compensation
% The so called phase should be in rad, the so-called angle should be in deg.
t=InputSignal.t;
% f=InputSignal.f;
y=InputSignal.y;

%% Directly perform the windowed FFT and Extract the amplitudes of harmonics
WindowTypes={'rect','hanning','blackmanharris'};
WindowFlag = find(strcmp(WindowName,WindowTypes));

if   isempty(WindowFlag) == 1
    error('We can not find this window in the prescribed window types!')
else
    switch WindowFlag
        case  1
            Wind=1;                         %% Time-domian window function
            Window_coeff=1;                 %% Amplitude correction coefficent to frequency-domain data
        case  2
            Wind=hanning(length(t))';       %% Time-domian window function
            Window_coeff=2;                 %% Amplitude correction coefficent to frequency-domain data
        case 3
            Wind=blackmanharris(length(t))';
            Window_coeff=2.7874566;         %% This coefficient is obtained using another code ; BlackmanHarrisAmplitudeCorrection.m
        otherwise
            error('Something wrong with the code')
    end
end

[Freq , Y ,AngY ] = MyFFTRobust(t , Wind.*y );   % PhaY in deg
% [~,     F ,AngF ] = MyFFTRobust(t , Wind.*f );
rawSpectrum.AngY=AngY;
rawSpectrum.Freq=Freq;


% % rawSpectrum.F=Window_coeff*F.*exp(1j*AngF*pi/180);
 rawSpectrum.Y=Window_coeff*Y.*exp(1j*AngY*pi/180);


