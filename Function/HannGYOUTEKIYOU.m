function [f,X] = HannGYOUTEKIYOU(t,x,fs)
% 步骤1：加窗
% 这里使用汉宁窗（Hanning Window），可以根据需要选择其他窗函数
x=x';
L = length(t);
N = 2^(nextpow2(L+1)-1); % 设置FFT点数
window = hann(N); % 创建汉宁窗，N 是信号的长度
x_windowed = x(1:N) .* window; % 将时域信号的前N个样本和窗函数相乘

% 步骤2：进行FFT变换
% 对加窗后的时域信号进行快速傅里叶变换得到频域信号
X = fft(x_windowed);

% 步骤3：计算频谱的幅值
X_mag = abs(X);

% 步骤4：修正幅值并进行对称性处理
X_mag_corrected = X_mag / (0.5 * N); % 对于汉宁窗，修正因子是1/（0.5*N）

% 对称性处理：将右半部分的幅值翻倍，合并到左半部分
X_mag_corrected(2:end) = 2 * X_mag_corrected(2:end);

% 计算频率轴，表示频谱中不同频率对应的频率值
% FFT的频率范围是从0到Fs，且频率间隔为 Fs/N
frequencies = (0:N-1) * fs / N;

% 只保留频谱的前一半（从0到Fs/2）
half_N = ceil(N/2);
frequencies = frequencies(1:half_N);
X_mag_corrected = X_mag_corrected(1:half_N);
X = X_mag_corrected;
f = frequencies;