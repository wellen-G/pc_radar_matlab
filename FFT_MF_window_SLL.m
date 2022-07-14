%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   匹配滤波器的FFT-IFFT的实现
%   设信号为时宽10us 带宽30MHz的LFM脉冲，采样率为80MHz。数学分析和仿真。
%   分析加窗后的旁瓣的变化 信噪比的损失
%   2021/12/1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
fs=8e7;               %采样频率为80MHz        采样频率为10MHz时，增益更接近理论值
Ts=1/fs;
B=30E6;              %LFM信号带宽为30MHz
T=10e-6;             %LFM信号时宽为10us
% N=2048;            %FFT的点数
SNR=-10;
k=B/T;               %调频斜率
A_xnoise=1;
n=round(T*fs);       %采样点个数
N=1024;                    %FFT的点数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
t=linspace(-T/2,T/2,n);
s_Preset=exp(1j*pi*k*t.^2);       %产生时宽为10us,带宽为30MHz的线性调频信号
xnoise_Preset= A_xnoise*randn(size(s_Preset))+1i*A_xnoise*randn(size(s_Preset));  %产生与信号等长度的噪声

if B<fs
    %%%%%     设计巴特沃斯低通滤波    %%%%%%%%%%%%%%%%%%%%%%%%%
    Wc=B/fs;                                          %截止频率为B/2Hz
    [b,a]=butter(20,Wc);
    xnoise=filter(b,a,xnoise_Preset);    %对噪声进行低通滤波
else
    xnoise=xnoise_Preset;
end
 %xnoise=xnoise_Preset;
%plot(t,s_Preset)

% s=s_Preset;SNR_in=10*log10(sum(s.*conj(s))/sum(xnoise.*conj(xnoise)));
% p_ss=10*log10(sum(s.*conj(s))/length(s));p_nn=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
[s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR,xnoise);   %得到给定信噪比下的信号以及噪声
h=conj(fliplr(s));      %匹配滤波起的冲激响应

[s_out,xnoise_out]=FFT_MF(s,xnoise,h);               %FFT实现匹配滤波 
[SNR_out,P_s_max,p_n_out]=SNR_max_out(s_out,xnoise_out);     %计算脉压后的峰值信噪比和信号峰值功率
Gain_MF=SNR_out-SNR_in;
%%%%%%%%%%%%%%%%%%%%% 考虑回波与发射信号不重合 %%%%%%%%%%%%%%%%%
N_return=40;
for i=1:N_return
    
    t_1=linspace(-T/2-Ts*(i/(N_return)),T/2-Ts*(i/(N_return)),n);
    s_reture_2(i,:)=exp(1j*pi*k*t_1.^2); 
    PS_reture_in_2(i)=10*log10(sum(s_reture_2(i,:).*conj(s_reture_2(i,:)))/length(s_reture_2(i,:)));
    SNR_in_reture_2(i)=PS_reture_in_2(i)-10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
%     [s_reture_2(i,:),s_add_noise_reture_2(i,:),SNR_in_reture_2(i),xnoise]=add_SNR_xnoise(s_return_2(i,:),SNR,xnoise);   %得到给定信噪比下的信号以及噪声
    [s_out_reture_2(i,:),xnoise_out]=FFT_MF(s_reture_2(i,:),xnoise,h);               %FFT实现匹配滤波 
    [SNR_out_reture_2(i),P_s_max_reture_2(i)]=SNR_max_out(s_out_reture_2(i,:),xnoise_out);     %计算脉压后的峰值信噪比和信号峰值功率
    Gain_MF_reture(i)=SNR_out_reture_2(i)-SNR_in_reture_2(i);
    
end

figure(6);   % 回波与发射信号不重合时
t_n=linspace(0,1,N_return);
plot(t_n,Gain_MF_reture,'linewidth',1),title('回波与发射信号不重合时,脉压增益与采样点的关系'),xlabel('相差采样点的个数'),ylabel('增益/dB');
clear i N_return;
    
t_1=linspace(-T/2-Ts/2,T/2-Ts/2,n);
s_return=exp(1j*pi*k*t_1.^2); 
[s_reture,s_add_noise_reture,SNR_in_reture,xnoise]=add_SNR_xnoise(s_return,SNR,xnoise);   %得到给定信噪比下的信号以及噪声
[s_out_reture,xnoise_out]=FFT_MF(s_reture,xnoise,h);               %FFT实现匹配滤波 
[SNR_out_reture,P_s_max_reture]=SNR_max_out(s_out_reture,xnoise_out);     %计算脉压后的峰值信噪比和信号峰值功率
Gain_MF_reture=SNR_out_reture-SNR_in_reture;


%%%%%%%%%%%%%%%%%%%%%% 进行加窗  %%%%%%%%%%%%%%%%%%%%%%%%%

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
s_out_addw=zeros(w_N,2^round(log2(2*length(s)-1)));       %脉压后的信号
xnoise_out_addw=zeros(w_N,2^round(log2(2*length(s)-1)));   %脉压后的噪声
SNR_out_addw=zeros(1,w_N);             %脉压后的峰值信噪比
P_s_max_addw=zeros(1,w_N);             %脉压后的信号峰值功率
SLL_addw=zeros(1,w_N);                 %脉压后的信号旁瓣电平
SLL_addw_interp=zeros(1,w_N);           %脉压后的信号旁瓣电平_做插值后
for i=1:w_N
  xnoise_addw(i,:)=xnoise.*(w(:,i)');
  s_addw(i,:)=s.*w(:,i)';           %加入汉宁窗
  h_addw=h.*w(:,i)';
  [s_out_addw(i,:),xnoise_out_addw(i,:)]=FFT_MF(s_addw(i,:),xnoise_addw(i,:),h);      %FFT实现匹配滤波 
  [SNR_out_addw(i),P_s_max_addw(i)]=SNR_max_out(s_out_addw(i,:),xnoise_out_addw(i,:));     % 得到输出信噪比以及信号峰值功率
  [SLL_addw(i)]=SLL_CAL(s_out_addw(i,:));     %计算旁瓣电平
  
end


%%%%%%%%%%%%%%%%%%   做插值运算 %%%%%%%%%%%%%%%%%%%%
N_interp=10;      %内插的数量   
s_out_addw_interp=zeros(w_N,2^round(log2(2*length(s)-1))*N_interp);       %脉压后的信号
for i=1:w_N
    s_out_addw_interp(i,:)=interp(s_out_addw(i,:),N_interp);
    [SLL_addw_interp(i)]=SLL_CAL(s_out_addw_interp(i,:));        %插值后的旁瓣电平
end
 
Ts_interp=Ts/N_interp;

[s_out_addw_peak_1,~,tmax_1,~]=Peak(s_out_addw_interp(1,:),Ts_interp);     %脉压后的信号极值点'
[s_out_addw_peak_2,~,tmax_2,~]=Peak(s_out_addw_interp(2,:),Ts_interp);     %脉压后的信号极值点'
[s_out_addw_peak_3,~,tmax_3,~]=Peak(s_out_addw_interp(3,:),Ts_interp);     %脉压后的信号极值点'
[s_out_addw_peak_4,~,tmax_4,~]=Peak(s_out_addw_interp(4,:),Ts_interp);     %脉压后的信号极值点'
[s_out_addw_peak_5,~,tmax_5,~]=Peak(s_out_addw_interp(5,:),Ts_interp);     %脉压后的信号极值点'
[s_out_addw_peak_6,~,tmax_6,~]=Peak(s_out_addw_interp(6,:),Ts_interp);     %脉压后的信号极值点'

%  
% figure(1)
% plot(t,real(h.*w(:,2)'))

%%%%%%%%%%%%%%%%%%%%%5 绘图  %%%%%%%%%%%%%%%%%%%%%%%%%%%%55
figure(1);   %  画出线性调频信号的时域和 频域波形
subplot(211);
plot(t,real(s_Preset)),title('LFM信号时域'),xlabel('t/s'),ylabel('幅度');
S=fftshift(fft(s_Preset,N));
f=linspace(-fs/2,fs/2,N);
subplot(212);
plot(f,abs(S)),title('LFM信号频谱'),xlabel('f/Hz'),ylabel('幅度');

figure(4);   %  画出线性调频信号的时域和 频域波形
subplot(211);
plot(t,real(xnoise)),title('噪声时域'),xlabel('t/s'),ylabel('幅度');
NK=fftshift(fft(xnoise,N));
f=linspace(-fs/2,fs/2,length(NK));
subplot(212);
plot(f,abs(NK)),title('噪声频谱'),xlabel('f/Hz'),ylabel('幅度');

figure(5);   % 回波与发射信号不重合时
hold on
plot(t,real(s_Preset),'-*','linewidth',1),title('LFM信号时域'),xlabel('t/s'),ylabel('幅度');
plot(t,real(s_return),'-+','linewidth',1),title('LFM信号时域'),xlabel('t/s'),ylabel('幅度'),ylim([-1.5,1.5]);
legend('发射信号','回波信号')
hold off


figure(3)
hold on
t_2=Ts:Ts:(length(s_out))*Ts;
plot(t_2,abs(s_out)),xlim([0.9E-5,1.1E-5]),title('匹配滤波后信号时域波形（FFT）')
plot(t_2,abs(xnoise_out)),xlim([0.9E-5,1.1E-5])
plot(t_2,abs(xnoise_out+s_out)),xlim([0.9E-5,1.1E-5])
hold off

figure(2)
t_interp=Ts_interp:Ts_interp:(length(s_out_addw_interp(1,:)))*Ts_interp;
m_1=0.95E-5;m_2=1.05E-5;
subplot(231)
hold on
plot(t_interp,20*log10(abs(s_out_addw_interp(1,:)))),xlim([m_1,m_2]),title('匹配滤波后信号时域波形（汉宁窗）')
plot(tmax_1,20*log10(abs(s_out_addw_peak_1)),'r+') 
hold off
subplot(232)
hold on
plot(t_interp,20*log10(abs(s_out_addw_interp(2,:)))),xlim([m_1,m_2]),title('匹配滤波后信号时域波形（切比雪夫窗）')
plot(tmax_2,20*log10(abs(s_out_addw_peak_2)),'r+')
hold off
subplot(233)
hold on
plot(t_interp,20*log10(abs(s_out_addw_interp(3,:)))),xlim([m_1,m_2]),title('匹配滤波后信号时域波形（三角窗）')
plot(tmax_3,20*log10(abs(s_out_addw_peak_3)),'r+')
hold off
subplot(234)
hold on
plot(t_interp,20*log10(abs(s_out_addw_interp(4,:)))),xlim([m_1,m_2]),title('匹配滤波后信号时域波形（海明窗）')
plot(tmax_4,20*log10(abs(s_out_addw_peak_4)),'r+')
hold off
subplot(235)
hold on
plot(t_interp,20*log10(abs(s_out_addw_interp(5,:)))),xlim([m_1,m_2]),title('匹配滤波后信号时域波形（布拉克曼窗）')
plot(tmax_5,20*log10(abs(s_out_addw_peak_5)),'r+')
hold off
subplot(236)
hold on
plot(t_interp,20*log10(abs(s_out_addw_interp(6,:)))),xlim([m_1,m_2]),title('匹配滤波后信号时域波形（未加窗）')
plot(tmax_6,20*log10(abs(s_out_addw_peak_6)),'r+')
hold off

 [SLL_addw_interp(5)]=SLL_CAL(s_out_addw_interp(5,:));  
 
function [s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR_Preset,xnoise)
    %  add_SNR_xnoise在信号中加入固定信噪比噪声
    %  s_Preset为待加入噪声的信号    SNR_Preset为信噪比  A_xnoise为噪声的标准差
    %  s为信号     s_add_noise 加入噪声后的信号   SNR_in输入信噪比    xnoise为噪声 
    Pn_add=sum(xnoise .*conj(xnoise));                % 计算加入噪声的能量 
    Ps_add=sum(s_Preset.*conj(s_Preset));             % 计算单位信号的能量 
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
   N3=length(xnoise_in);
   N_s=2^ceil(log2(N1+N2-1));
   SK=fft(s_in,N_s);
   HK=fft(h_in,N_s);
   SK_out=SK.*HK;          %LFM信号与MF做FFT然后相乘
   s_mf_out=ifft(SK_out,N_s);
  
   N_n=2^ceil(log2(N3+N2-1));
   NK=fft(xnoise_in,N_n);
   HK_n=fft(h_in,N_n);
   NK_out=NK.*HK_n; 
   xnoise_mf_out_temp=ifft(NK_out,N_n);
   xnoise_mf_out=xnoise_mf_out_temp(N_n/2-N_n/2*(N_s/N_n)+1:N_n/2+N_n/2*(N_s/N_n));
   
end

function  [SNR_out,p_s_max,p_n]=SNR_max_out(s,xnoise)
  %求输出信噪比及峰值功率
  p_s_max=10*log10(max(s)*conj(max(s)));
  p_n=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
  SNR_out =p_s_max- p_n;
end


function  [vmax,vmin,tmax,tmin]=Peak(signal,Ts)
   %  求序列的极值
    t=Ts:Ts:(length(signal))*Ts;
    Lmax=diff(sign(diff(abs(signal))))==-2;
    Lmin=diff(sign(diff(abs(signal))))==2;
    Lmax=[false,Lmax,false];
    Lmin=[false,Lmin,false];
    tmax=t(Lmax);
    tmin=t(Lmin);
    vmax=signal(Lmax);
    vmin=signal(Lmin);
end

function [SLL]=SLL_CAL(signal)
    %求信号的旁瓣电平
   
%     signal=signal(round(11/30*length(signal)):round(19/30*length(signal)));
    Lmax=diff(sign(diff(abs(signal))))==-2;    %得到极值点
    Lmax=[false,Lmax,false];
    vmax=signal(Lmax);
    [S_max,S_max_index]=max(vmax);             %求出主瓣电平和位置
    S_sec_sum=[vmax(S_max_index-5:S_max_index-1),vmax(S_max_index+1:S_max_index+5)];  %得到主瓣周围的旁瓣
    S_sec=max(S_sec_sum);               %得到最大旁瓣
%     Signal_sort=sort(abs(vmax));    %对序列排序
%     S_max= Signal_sort(length(Signal_sort));
%     S_sec=Signal_sort(length(Signal_sort)-1);
    SLL=20*log10(abs(S_sec)/abs(S_max));
 end
