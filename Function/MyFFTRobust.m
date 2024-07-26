%--------------------------------------%
%   对给定信号进行FFT分析,并绘制图形     %
%   参数: 时间向量t, 信号向量x          %
%   要求: t等间隔分布, x与t同维         %
%   Programmed by Xing Wang              %
%   2019-10-31                        %
%--------------------------------------%
function [f, Xabs, Pha] = MyFFTRobust(t, x)
size_t = size(t);
if(size_t(1) > size_t(2))
    t = t';
end
size_x = size(x);
if(size_x(1) > size_x(2))
    x = x';
end
if t(1)~=0
    error('Please shift the time to start from 0!')
end


FS = length(t) / t(end);
N = length(x);
X = fft(x);


if mod(N,2) == 0 % number is even 偶数
 f = linspace(0,FS/2,N/2+1);
 Xabs = 2*abs(X(1:N/2+1)) /N;
 Pha = (angle(X(1:N/2+1)))*180/pi;
 
else %number is odd   
f = [0:(N-1)/2]*FS/N;
Xabs = 2*abs(X(1:((N-1)/2+1)) )/N;
Pha = (angle(X(1:((N-1)/2+1))))*180/pi;
end

% figure
% subplot(2, 1, 1);
% plot(t, x, '.-');
% xlabel('t(s)');
% ylabel('x(t)');
% grid on;
% subplot(2, 1, 2);
% semilogx(f, Xabs, '.-');
% xlabel('Frequency(Hz)');
% ylabel('|X(f)|');
% grid on;
