%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  设噪声带宽是5/10MHz，单点频复信号频率是500KHz，采样率是5MHz，仿真分析FFT的增益。
%   SNR可调   FFT点数可调
% 2021/11/30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc
N=1024;        % FFT的点数为   N=1024
fs=10E6;        % 采样率       fs=5MHz   
Ts=1/fs;
%     此时分辨率为fs/N= 4.8828e+03    4.8KHz
B_noise=10e6;   %噪声带宽
fd=1.25e6;     % 信号频率     fd=1.25MHz
% fd=1.25e6+D(i_avr)*fs/N;     % 信号频率     fd=1.25MHz
% C=0.5;
% fd=1.25e6+C*fs/N;
SNR_Preset=-10;         %初始信噪比，以dB为单位
A_xnoise=1;                  %噪声标准差
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:Ts:N*Ts-Ts;
n=length(t); %采样点的个数
s_Preset=exp(1i*(2*pi*fd*t+pi/6));      %先定义幅值为1单频带信号
% s_Preset=exp(1i*(2*pi*fd*t+pi/6))+exp(1i*(2*pi*1.25e6*t+pi/6)); 
xnoise_Preset= A_xnoise*randn(size(s_Preset))+1i*A_xnoise*randn(size(s_Preset));  %产生与信号等长度的噪声
if B_noise<fs
    %%%%%     设计巴特沃斯低通滤波    %%%%%%%%%%%%%%%%%%%%%%%%%
    Wc=B_noise/fs;                                          %截止频率为B/2Hz
    [b,a]=butter(20,Wc);
    xnoise=filter(b,a,xnoise_Preset);    %对噪声进行低通滤波
else
    xnoise=xnoise_Preset;
end

%%%%%%%%%%%%%%  信号不变时噪声带宽的变化对输入信噪比的影响  %%%%%%%%%%%%%%55
[SNR_in,p_s,p_n]=SNR_caj_time(s_Preset,xnoise);   
% s=s_Preset;
[s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR_Preset,xnoise);
%%%%%%%%%%%%%%%%%%%%  验证能量守恒  %%%%%%%%%%%%%%%%%%%%%%%
[SNR_in,p_s,p_n]=SNR_caj_time(s,xnoise);   
SK=fftshift(fft(s,N));
NK=fftshift(fft(xnoise,N)); 
% XK=fftshift(fft(s_add_noise,N)); 
[SNR_K_out,P_SK,P_NK]=SNR_caj_fre(SK,NK); 

[SNR_FFT_out,p_s_max,p_NK_2]=SNR_max_out(SK,NK,fs,B_noise);  % FFT后的峰值信噪比
Gain_FFT=SNR_FFT_out-SNR_in;
%%%%%%%%%%%%%%%%%% 频点位置和FFT后输出信噪比的关系 %%%%%%%%%%%%%%%%%%%%%%%%%%
%  调用avr_SNR程序来执行

%%%%%%%%%%%%%%%%%%%%% 不同窗函数对输出信噪比的影响  %%%%%%%%%%%%%%%%%   


 w(:,1)=hann(length(t));      %汉宁窗函数
 w(:,2)=chebwin(length(t));   %切比雪夫窗函数
%   wvtool(w(:,2))
 w(:,3)=triang(length(t));    %三角窗函数
 w(:,4)=hamming(length(t));   %海明窗函数
 w(:,5)=blackman(length(t));  %布拉克曼窗函数
 w(:,6)=boxcar(length(t));     %矩形窗函数（相当于未加窗)
  
w_N=length(w(1,:));                    %加入窗函数的数量
xnoise_addw=zeros(w_N,length(s));      %加窗后的噪声
s_addw=zeros(w_N,length(s));           %加窗后的信号
SK_addw=zeros(w_N,N);                  %FFT后的信号
NK_addw=zeros(w_N,N);                  %FFT后的噪声
SNR_in_addw=zeros(1,w_N);              %FFT前的峰值信噪比
SNR_out_addw=zeros(1,w_N);             %FFT后的峰值信噪比
Gain_in_addw=zeros(1,w_N);             %FFT后的峰值信噪比增益
P_s_max_addw=zeros(1,w_N);             %FFT后的信号峰值功率

for i=1:w_N
  xnoise_addw(i,:)=xnoise.*(w(:,i)');    %噪声加窗
  s_addw(i,:)=s.*w(:,i)';                %信号加窗
  [SNR_in_addw(1,i),~,~]=SNR_caj_time(s_addw(i,:),xnoise_addw(i,:));   
  SK_addw(i,:)=fftshift(fft(s_addw(i,:),N));
  NK_addw(i,:)=fftshift(fft(xnoise_addw(i,:),N)); 
  [SNR_out_addw(1,i),~,~]=SNR_max_out( SK_addw(i,:),NK_addw(i,:),fs,B_noise);  % FFT后的峰值信噪比
  Gain_in_addw(1,i)=SNR_out_addw(1,i)-SNR_in_addw(1,i);          %FFT后的峰值信噪比增益
end


%%%%%%%%%%%%%%  绘图 %%%%%%%%%%%%%%%%%%%%%%%
% f=linspace(-fs/2,fs/2-fs/length(SK),length(SK));
% figure(1)   %%%%%%%%%  信号  %%%%%%%%%%
% % subplot(311)
% % plot(t,real(s)), xlabel('t/s'),ylabel('幅度'),title('信号时域波形'),grid on;
% 
% subplot(211)
% plot(t,real(s),'linewidth',1), xlabel('t/s'),ylabel('幅度'), xlim([0 ,200*Ts]),title('信号时域波形'),grid on;
% 
% subplot(212)
% plot(f,abs(SK),'linewidth',1), xlabel('f/Hz'),ylabel('幅度'),title('信号频域域波形'),grid on;
% 
% % figure(2)        %%%%%%%%%%%  噪声  %%%%%%%%%%
% % subplot(211)
% % plot(t,abs(xnoise)), xlabel('t/s'),ylabel('幅度'),title('噪声时域波形'),grid on;
% % 
% % subplot(212)
% % 
% % plot(f,abs(NK),'linewidth',1), xlabel('f/Hz'),ylabel('幅度'),title('噪声频域域波形'),grid on;
% 
% figure(3)        %%%%%%%%%%%  信号+噪声  %%%%%%%%%%
% % subplot(311)
% % plot(t,real(s_add_noise)), xlabel('t/s'),ylabel('幅度'),title('信号+噪声时域波形'),grid on;
% subplot(211)
% plot(t,real(s_add_noise),'linewidth',1), xlabel('t/s'),ylabel('幅度'),xlim([0 ,200*Ts]),title('信号+噪声时域波形'),grid on;
% subplot(212)
% plot(f,abs(XK),'linewidth',1), xlabel('f/Hz'),ylabel('幅度'),title('信号+噪声频域域波形'),grid on;
% 
figure(4)  
subplot(211)
plot(f,abs(SK),'linewidth',1.2), xlabel('f/Hz'),ylabel('幅度'),title('信号频域域波形'),grid on;
subplot(212)
stem(f,abs(SK) ,'linewidth',1.2), xlabel('f/Hz'),ylabel('幅度'),xlim([1.2e6 ,1.3e6]),title('信号频域域波形(小区间内)'),grid on;



function [s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR_Preset,xnoise)
    %  add_SNR_xnoise在信号中加入固定信噪比噪声
    %  s_Preset为待加入噪声的信号    SNR_Preset为信噪比  A_xnoise为噪声的标准差
    %  s为信号     s_add_noise 加入噪声后的信号   SNR_in输入信噪比    xnoise为噪声 
    Pn_add=sum(xnoise .*conj(xnoise ));                  % 计算加入噪声的功率 
    Ps_add=sum(s_Preset.*conj(s_Preset));                % 计算单位信号的功率 
    A_1= sqrt((10^(SNR_Preset/10))*(Pn_add/Ps_add));     %计算对应信噪比下信号的幅度
    s=A_1.*s_Preset;                                     %得到对应信噪比下的信号
    s_add_noise=s+xnoise;                                % 加入噪声的信号
    Ps_add_out=sum(s.*conj(s));
    SNR_in=10*log10(Ps_add_out/Pn_add);                   %输入信噪比
end
function  [SNR_out,p_s,p_n]=SNR_caj_fre(s,xnoise)
  %求信噪比及平均功率
  p_s=10*log10(sum(s.*conj(s))/length(s));
  p_n=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
  SNR_out =p_s- p_n;
end
function  [SNR_out,p_s,p_n]=SNR_caj_time(s,xnoise)
  %求信噪比及能量
  p_s=10*log10(sum(s.*conj(s)));
  p_n=10*log10(sum(xnoise.*conj(xnoise)));
  SNR_out =p_s- p_n;
end
function  [SNR_out,p_s_max,p_n]=SNR_max_out(s,xnoise,fs,B_n)
  %求输出信噪比及峰值功率
  %fs为采样率  B_n为噪声带宽
  n=length(xnoise);
  n1=round(n/2-n/2*(B_n/fs))+1;
  n2=round(n/2+n/2*(B_n/fs));
  p_s_max=10*log10(max(s)*conj(max(s)));
  p_n=10*log10(sum(xnoise(n1:n2).*conj(xnoise(n1:n2)))/(n2-n1+1));
  SNR_out =p_s_max- p_n;
end