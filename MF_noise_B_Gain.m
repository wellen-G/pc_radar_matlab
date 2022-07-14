%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%理论分析和matlab仿真等功率谱密度的高斯白噪声，反映不同带宽下的噪声功率变化；
%设信号为时宽10us带宽30MHz/15MHz的LFM脉冲，采样率为80MHz；反映信号带宽对脉压增益的影响。
% 2021/11/30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
fs=8e7;                    %采样频率为80MHz
Ts=1/fs;
B=15E6;                    %LFM信号带宽为30/15MHz
T=10e-6;                   %LFM信号时宽为10us
N=1024;                    %FFT的点数
SNR_Preset=-10;              %输入信号信噪比
k=B/T;                     %调频斜率
A_xnoise=1;                  %噪声标准差
n=round(T*fs);             %采样点个数
%
%    理想的脉压增益：BT      脉压的增益与输入信号的时宽、带宽有关
%       Gain=10*log10(BT)
%       10*log10(300)=24.7712     10*log10(150)=21.7609
%    而改变信号的带宽仅能改变输入信噪比，改变不了输出信噪比
%    脉压的输出信噪比与时宽有关
%    
%
%    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=linspace(-T/2,T/2,n);
s_Preset=exp(1j*pi*k*t.^2);          %产生时宽为10us,带宽为30MHz的线性调频信号
h=conj(fliplr(s_Preset));            %匹配滤波器的冲激响应
xnoise_Preset= A_xnoise*randn(size(s_Preset))+1i*A_xnoise*randn(size(s_Preset));  %产生与信号等长度的噪声

if B<fs
    %%%%%     设计巴特沃斯低通滤波    %%%%%%%%%%%%%%%%%%%%%%%%%
    Wc=B/fs;                                          %截止频率为B/2Hz
    [b,a]=butter(20,Wc);
    xnoise=filter(b,a,xnoise_Preset);    %对噪声进行低通滤波
else
    xnoise=xnoise_Preset;
end


 
%%%%%%%%%%%%%%%%%%%%%%%%  信号带宽与输出信噪比的影响  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%    噪声功率谱密度不变  %%%%%%%%%%%%%%%%%%%%%%   
[SNR_in,p_s,p_n]=SNR_caj(s_Preset,xnoise);
[s_mf_out,xnoise_mf_out]=FFT_MF(s_Preset,xnoise,h);
[SNR_out,p_s_max]=SNR_max_out(s_mf_out,xnoise_mf_out);
%%%%%%%%%%%%%%%%%%%%    固定输入信噪比  %%%%%%%%%%%%%%
[s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR_Preset,xnoise);
h=conj(fliplr(s));  %匹配滤波起的系统函数
[SNR_in_2,p_s,p_n]=SNR_caj(s,xnoise);
[s_mf_out,xnoise_mf_out]=FFT_MF(s,xnoise,h);
[SNR_out,p_s_max]=SNR_max_out(s_mf_out,xnoise_mf_out);

%%%%%%%%%%%%%%% 绘图 %%%%%%%%%%%%%%%%%%%%%%%%%%%%5 
figure(1);   %  画出线性调频信号的时域和 频域波形
subplot(211);
plot(t,real(s_Preset)),title('LFM信号时域'),xlabel('t/s'),ylabel('幅度');
S=fftshift(fft(s_Preset,N));
f=linspace(-fs/2,fs/2,N);
subplot(212);
plot(f,abs(S)),title('LFM信号频谱'),xlabel('f/Hz'),ylabel('幅度');

figure(2);   %  画出线性调频信号的时域和 频域波形
subplot(211);
plot(t,real(xnoise)),title('噪声时域'),xlabel('t/s'),ylabel('幅度');
NK=fftshift(fft(xnoise,N));
f=linspace(-fs/2,fs/2,length(NK));
subplot(212);
plot(f,abs(NK)),title('噪声频谱'),xlabel('f/Hz'),ylabel('幅度');




function [s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR_Preset,xnoise)
    %  add_SNR_xnoise在信号中加入固定信噪比噪声
    %  s_Preset为待加入噪声的信号    SNR_Preset为信噪比  A_xnoise为噪声的标准差
    %  s为信号     s_add_noise 加入噪声后的信号   SNR_in输入信噪比    xnoise为噪声
    Pn_add=sum(xnoise .*conj(xnoise ));                % 计算加入噪声的能量 
    Ps_add=sum(s_Preset.*conj(s_Preset));               % 计算单位信号的能量 
    A_1= sqrt((10^(SNR_Preset/10))*(Pn_add/Ps_add));     %计算对应信噪比下信号的幅度
    s=A_1.*s_Preset;                                    %得到对应信噪比下的信号
    s_add_noise=s+xnoise;                            % 加入噪声的信号
    Ps_add_out=sum(s.*conj(s));
    SNR_in=10*log10(Ps_add_out/Pn_add);                   %输入信噪比
end


function   [s_mf_out,xnoise_mf_out]=FFT_MF(s_in,xnoise_in,h_in)
   % FFT实现匹配滤波
   
   xnoise_in=[xnoise_in,xnoise_in,xnoise_in,xnoise_in];   %噪声序列延长
  
   N1=length(s_in);
   N2=length(h_in);
   N_sum=N1+N2-1;
   v=ceil(log2(N_sum));
   N_s=2^v;
   SK=fft(s_in,N_s);
   HK=fft(h_in,N_s);
   SK_out=SK.*HK;          %LFM信号与MF做FFT然后相乘
   s_mf_out=ifft(SK_out,N_s);
   
   N3=length(xnoise_in);
   N_n=2^ceil(log2(N1+N3-1));
   NK=fft(xnoise_in,N_n);
   HK_n=fft(h_in,N_n);
   NK_out=NK.*HK_n;
   xnoise_mf_out_temp=ifft(NK_out,N_n);
   xnoise_mf_out=xnoise_mf_out_temp(N_n/2-N_n/2*(N_s/N_n)+1:N_n/2+N_n/2*(N_s/N_n));
   
   
end

function  [SNR_out,p_s_max]=SNR_max_out(s,xnoise)
  %求输出信噪比及峰值功率
  p_s_max=10*log10(max(s)*conj(max(s)));
  p_n=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
  SNR_out =p_s_max- p_n;
end

function  [SNR_out,p_s,p_n]=SNR_caj(s,xnoise)
  %求信噪比及平均功率
  p_s=10*log10(sum(s.*conj(s))/length(s));
  p_n=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
  SNR_out =p_s- p_n;
end




