clc 
clear 
close all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_length=127;                     %码长为127
fm=72e6;                          %码频为72MHz
fc=10e9;                          %信号载频为 10GHz
fs=72e6*4;                         %采样率为200MHz
Ts=1/fs;
Duty_Ratio= 0.10;                  %占空比为10%
PRT=m_length*(1/fm)/Duty_Ratio;    %脉冲重复周期
PRF=1/PRT;                        %脉冲重复频率
n_add=round(100e-3/(m_length/fm)*Duty_Ratio);        %相参累计个数
n_add=2^round(log2(n_add));

 n_add=256;
CPI=n_add*PRT;                       %相参处理时间
SNR_Preset=0;                       %输入信噪比
A_xnoise =1;                            %噪声标准差
N_PRT=m_length*(fs/fm)/Duty_Ratio;    %一个周期内的采样点数
T_width=m_length/fm;              %信号的时宽
B=fm;
%%%%%%%%%%%%%%%%%%%%%%  产生m序列   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 高斯随机信号 
u = idinput(m_length, 'prbs',[0,1],[-1,1])';   %产生m序列
%  figure(12) ,stairs(u,'linewidth',1) ,title('m序列')
[c,lags] = xcorr(u);    % 归一化后m序列的自相关函数
%  figure(11),  plot(lags,c,'linewidth',1),title('m序列自相关')   %绘制m序列的自相关
 m_sum=[];
 z=zeros(m_length*10*(1-Duty_Ratio),1);
 for i=1:n_add  
    m_sum=[m_sum,idinput(m_length, 'prbs',[0,1],[-1,1])'];%,z'];   
 end
[c_sum,lags_sum] = xcorr(u ,m_sum);    % 归一化后m序列的循环自相关函数
 % figure(10) , plot(c_sum,'linewidth',1),title('m序列的循环自相关')     %绘制m序列的循环自相关
%%%%%%%%%%%%%%%%%%%%%%%% 产生伪随机相位编码调制信号 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

Phase_code_signal=Phase_code_modulation( u,fm,fc,fs,Duty_Ratio,n_add) ;      %产生发射信号
%    plot( real(Phase_code_signal)),xlim([0,2e5]),ylim([-2,2]),ylim([-2,2]), title('伪随机相位编码信号') 
s_length=length(Phase_code_signal);                        %发射信号的长度
xnoise_Preset= A_xnoise*randn(size(Phase_code_signal))+1i*A_xnoise*randn(size(Phase_code_signal));  
 [signal_demodulation]=demodulation(Phase_code_signal,fc,fs);     %对发射信号解调 
  h=conj(fliplr(signal_demodulation(1:m_length*round(fs/fm))));   %产生匹配滤波器而对冲激响应 
 

 SK_mf=fftshift(fft(signal_demodulation));
 f_s=linspace(-fs/2,fs/2,length(Phase_code_signal));
 % figure(1),plot(f_s,abs(SK_mf));

 %%%%%%%%% 验证脉压的增益  %%%%%%%%%%%%%%%%%%%%%%%
   %%%%%     设计巴特沃斯低通滤波    %%%%%%
  Wc=B/fs;                                          %截止频率为B/2Hz
  N_filter=20;     % 滤波器的阶数
  [b_filter,a_filter]=butter(N_filter,Wc);

  xnoise=filter(b_filter,a_filter,xnoise_Preset(1:N_PRT/10));    %对噪声进行低通滤波
%   s=filter(b_filter,a_filter,signal_demodulation(1:N_PRT/10));     %对信号进行低通滤波
  [ signal_filter_return]=filter_compose(signal_demodulation,a_filter,b_filter,N_PRT,Duty_Ratio,n_add);
  [ xnoise_filter_return]=filter_compose(xnoise_Preset,a_filter,b_filter,N_PRT,Duty_Ratio,n_add);
  
  %%%%%%%%%%%%%%  分析脉压的增益  %%%%%%%%%%%%%%%
  %  以第一个周期为例  ，无延时   无多普勒效应
  [s_snr,s_add_noise,SNR_in,xnoise_snr]=add_SNR( signal_filter_return(1:N_PRT/10),SNR_Preset,xnoise_filter_return(1:N_PRT/10));
  
  [signal_mf_out,noise_mf_out]=FFT_MF(s_snr(1:N_PRT/10),xnoise_snr(1:N_PRT/10),h);
  [SNR_in_MF,~,~]=SNR_caj_time(s_snr,  xnoise_snr);
  [SNR_out_MF,~,~]=SNR_max_out(signal_mf_out,noise_mf_out);
  Gain_MF=SNR_out_MF-SNR_in_MF;
 
  NK_mf=fftshift(fft(xnoise));
  f_n=linspace(-fs/2,fs/2,length(xnoise));
 % figure(2),plot(f_n,abs(NK_mf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%  产生回波信号 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%速度分辨率  
v_distinguish=3e8*PRF/2/n_add/fc;
%距离分辨率  3e8/fm/2
R_distinguish=3e8/fm/2;
R_max=(3e8/2)*PRT;     %最大不模糊距离
v_max=3e8*PRF/2/fc;    %最大不模糊速度

R_target=[1000,1000];               %目标的距离    R_max=(3e8/2)*PRT*(1-Duty_Ratio)=2381 R_min=(3e8/2)*PRT*Duty_Ratio=264
v_target=[400,435];                          %目标的速度    v_max=3e8*PRF/2/fc=850
A_target=[1,0.04];                        %目标回波信号的幅值衰减

%目标一的雷达回波
Phase_code_signal_return_1=Phase_code_signal_return_modulation(Phase_code_signal,fs,fc,Duty_Ratio,R_target(1),v_target(1),A_target(1),n_add);
%目标二的雷达回波
Phase_code_signal_return_2=Phase_code_signal_return_modulation(Phase_code_signal,fs,fc,Duty_Ratio,R_target(2),v_target(2),A_target(2),n_add);

%对回波信号解调
 [Phase_code_signal_return_1]=demodulation(Phase_code_signal_return_1,fc,fs);
 [Phase_code_signal_return_2]=demodulation(Phase_code_signal_return_2,fc,fs);
 xnoise=filter(b_filter,a_filter,xnoise_Preset);    %对噪声进行低通滤波
 % plot( real(Phase_code_signal_return_1)),xlim([0,2e5]),ylim([-2,2]), title('无多普勒效应时的回波信号')
 
 % plot(abs(fft(Phase_code_signal_return_2)))
%加入固定信噪比后的回波信号
%合成回波+噪声
[Phase_code_signal_return_sum,Phase_code_signal_return_sum_add_noise,SNR_in_sum_return,xnoise_sum_return]=add_SNR_xnoise_twotarget(Phase_code_signal_return_1,Phase_code_signal_return_2,SNR_Preset,A_xnoise,m_length*round(fs/fm));
%目标一回波+噪声
[Phase_code_signal_return_1 ,Phase_code_signal_return_1_noise, SNR_in_return(1),xnoise_return_1]=add_SNR_xnoise(Phase_code_signal_return_1,SNR_Preset,xnoise,m_length*round(fs/fm),N_filter);
%目标二回波+噪声
% [Phase_code_signal_return_2 ,Phase_code_signal_return_2_noise, SNR_in_return(2),xnoise_return_2]=add_SNR_xnoise(Phase_code_signal_return_2,SNR_Preset,A_xnoise,m_length*round(fs/fm));
 
% plot( real(Phase_code_signal_return_1_noise)), title('加入噪声后的回波信号'),xlim([0,2e5]),ylim([-2,2])
% plot( real(Phase_code_signal_return_1)),xlim([0,2e5]),ylim([-2,2]), title('无多普勒效应时的回波信号')

 SK_in=fftshift(fft(Phase_code_signal_return_1));
 SK_MF_in=fftshift(fft(xnoise_return_1));
 %脉压前后的带宽
%  figure(8), subplot(211), plot(linspace(-fs/2,fs/2,length(SK_in)),abs(SK_in)),title('信号频谱'),subplot(212), plot(linspace(-fs/2,fs/2,length( SK_MF_in)), abs(SK_MF_in)),title('噪声频谱');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   进行脉压 分析脉压  时宽  带宽  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [Phase_code_signal_return_1_mf_out,Phase_code_noise_return_1_mf_out]=FFT_MF(Phase_code_signal_return_1, xnoise_return_1,h);
%  [Phase_code_signal_return_2_mf_out,Phase_code_noise_return_2_mf_out]=FFT_MF(Phase_code_signal_return_2,xnoise_return_2,h);
 [Phase_code_signal_return_sum_mf_out,Phase_code_noise_return_sum_mf_out]=FFT_MF(Phase_code_signal_return_sum,xnoise_sum_return,h);
 
  % plot(real(Phase_code_signal_return_1_mf_out)), title('脉压后的信号'),xlim([0,2e5]),ylim([-2,2])
 
 %////////////////////  脉压后的时宽    取前五个脉冲时宽分析  /////////////////////////////
 signal_mf_width_analysis=Phase_code_signal_return_1_mf_out(1:N_PRT*5);
 place_max=zeros(1,5);  %最大值所在位置
 place_4dB=zeros(1,5);  %下降4dB所在位置

 N_interp=30;
%  s_interp = interp1( 1:1:N_PRT*5,signal_mf_width_analysis, 1/N_interp:1/N_interp:N_PRT*5,'spline')/max(signal_mf_width_analysis);
s_interp=interp(signal_mf_width_analysis, N_interp)/max(signal_mf_width_analysis);  %进行插值并归一化
  for i=1:5
     [~,place_max(i)]=max(s_interp((i-1)*N_PRT* N_interp+1:i*N_PRT* N_interp));
     place_max(i)=place_max(i)+(i-1)*N_PRT *N_interp;
 end
  for i=1:5   %找出下降4dB所在位置
     [~,place_4dB(i)]=min(abs(10*log10(abs(s_interp((i-1)*N_PRT* N_interp+1:i*N_PRT* N_interp)))+4));
     place_4dB(i)=place_4dB(i)+(i-1)*N_PRT* N_interp;
  end
 T_mf=sum(abs(place_max- place_4dB))*2/fs/N_interp/5;   %脉压后的时宽
 D_MF=T_width/T_mf;                                    %压缩比
 
 t_interp=0:1/fs/N_interp:(length(s_interp)-1)/fs/N_interp;
%脉压后的时宽
%  figure(9), plot(t_interp, 10*log10(abs(s_interp))),  xlim([(place_max(1)-150)*1/fs/N_interp,(place_max(1)+150)*1/fs/N_interp])
%//////////////////////// 脉压后的带宽/////////////////////////

 SK=fftshift(fft(Phase_code_signal_return_1));
 SK_MF=fftshift(fft( Phase_code_signal_return_1_mf_out));
 %脉压前后的带宽
%  figure(8), subplot(211), plot(linspace(-fs/2,fs/2,length(SK)),abs(SK)),title('脉压前频谱'),subplot(212), plot(linspace(-fs/2,fs/2,length( SK_MF)), abs(SK_MF)),title('脉压后频谱');
%  经过脉压后带宽不变为码频fm
 %%%%%%%%%%%%%%%%%%%%%%%%%%%  距离门重排  %%%%%%%%%%%%%%%%%%%%%%%%%%5 
 n_y=n_add;                                %相参累计个数
 n_x=m_length*(fs/fm)/Duty_Ratio;          %单个周期采样个数
 signal_return_1_array=zeros(n_x,n_y);
 noise_return_1_array=zeros(n_x,n_y);
%  signal_return_2_array=zeros(n_y,n_x);
 signal_return_sum_array=zeros(n_x,n_y);

for i=1:n_y
    for j=1:n_x
        signal_return_1_array(j,i)=Phase_code_signal_return_1_mf_out((i-1)*n_x+j);
        noise_return_1_array(j,i)=Phase_code_noise_return_1_mf_out((i-1)*n_x+j);
        signal_return_sum_array(j,i)=Phase_code_signal_return_sum_mf_out((i-1)*n_x+j);
    end
end
% %/////////////   距离 分辨率 ////////////////////////
% hold on
% for i=1:1:n_y
%     
%     plot(abs(signal_return_sum_array(:,i)))
% end
% hold off

% figure(16)  %距离门上叠加波形
% hold on
% for i=1:1:n_y
%     
%     plot(abs(signal_return_1_array(:,i)+noise_return_1_array(:,i)))
% end
% hold off

% figure(7);mesh(1:n_y,1:n_x,real(signal_return_1_array));
% figure(6);mesh(1:n_y,1:n_x,abs(signal_return_1_array+  noise_return_1_array));
% figure(5);mesh(1:n_y,1:n_x,real(signal_return_sum_array));

%%%%%%%%%%%%%%%%%%%%%%%  做FFT  得到速度分量%%%%%%%%%%%%%%%%%%%%%
signal_return_1_array_FFT=zeros(n_x,n_y);
noise_return_1_array_FFT=zeros(n_x,n_y);
% signal_return_2_array_FFT=zeros(n_x,n_y);
signal_return_sum_array_FFT=zeros(n_x,n_y);
for i=1:n_x
    signal_return_1_array_FFT(i,:)=fft(signal_return_1_array(i,:),n_add);
    noise_return_1_array_FFT(i,:)=fft(noise_return_1_array(i,:),n_add);
    signal_return_sum_array_FFT(i,:)=fft(signal_return_sum_array(i,:),n_add);
end

%%%%%%%%%%%%%%%%%%%%%%%%% 分析FFT的增益  %%%%%%%%%%%%%%%
[~, b_n]=max(signal_return_1_array(:,1));  %得到脉压输出最大点的位置
s_beforeFTT=signal_return_1_array(b_n,:);
n_beforeFTT=noise_return_1_array(b_n,:);
s_FFT_max=signal_return_1_array_FFT(b_n,:);  %得到信号FFT输出最大点的行
n_FFT_max=noise_return_1_array_FFT(b_n,:);   %得到噪声FFT输出最大点的行

% [~, n_fft]=max(s_FFT_max);

%%%%fft后的带宽
signal_return_1_array_FFT_interp=interp(signal_return_1_array_FFT(b_n,:), 50);
f_fft=linspace(0,fs,length(signal_return_1_array_FFT_interp));
f_fft_2=linspace(0,PRF,length(signal_return_1_array_FFT_interp));
% plot(f_fft,10*log10(abs(signal_return_1_array_FFT_interp)),'linewidth',1);
% plot(f_fft_2,10*log10(abs(signal_return_1_array_FFT_interp)),'linewidth',1);

% 理论上的增益为10*log10(N)  N为FFT 的点数  这里为脉冲的个数 n_add
[SNR_FFT_IN,~,~]=SNR_caj_time(s_beforeFTT,n_beforeFTT);
[SNR_FFT_out,~,~]=SNR_max_out(s_FFT_max,n_FFT_max);
Gain_fft=SNR_FFT_out-SNR_FFT_IN;

%%%%求距离门和速度门的位置
sort_A = sort(signal_return_sum_array(:,1), 'descend'); %降序排列
R_first_max = sort_A(1);
R_second_max = sort_A(2);
[~, R1_row] = find(signal_return_sum_array(:,1).' == R_first_max);    %距离门
[~, R2_row] = find(signal_return_sum_array(:,1).' == R_second_max);   

sort_B = sort(signal_return_sum_array_FFT(R1_row,:), 'descend'); %降序排列
v_first_max = sort_B(1);
v_second_max = sort_B(2);
[~, v1_row] = find(signal_return_sum_array_FFT(R1_row,:) == v_first_max);     %速度门
[~, v2_row] = find(signal_return_sum_array_FFT(R1_row,:) == v_second_max );    

% %/////////////   距离 分辨率 ////////////////////////
% figure(17)
% hold on
% for i=1:1:n_y
%     
%     plot(abs(signal_return_sum_array_FFT(:,i))),xlim([2000,3000])
% end
% hold off

%/////////////  速度分辨率 ////////////////////////
figure(4)
hold on
for i=1:1:n_x
    
    plot(abs(signal_return_sum_array_FFT(i,:))),xlim([100,140])
end
hold off




% figure(15)  %速度门上叠加波形
% hold on
% for i=1:1:n_x
%     
%     plot(abs(signal_return_1_array_FFT(i,:)+noise_return_1_array_FFT(i,:)))
% end
% hold off



figure(1);mesh(1:n_y,1:n_x,abs(signal_return_1_array_FFT));
figure(2);mesh(1:n_y,1:n_x,real(signal_return_1_array_FFT+ noise_return_1_array_FFT));
figure(3);mesh(1:n_y,1:n_x,abs(signal_return_sum_array_FFT));



function  [Phase_code_signal_return]=Phase_code_signal_return_modulation(Phase_code_signal,fs,fc,Duty_Ratio,R_target,v_target,A_target,n_add)
    %产生雷达回波信号
    %对发射信号延时，加多普勒频率
    c=3e8;    %光速
    n_delay=round(2*R_target*fs/c);    %由于距离引起的延时
    Phase_code_signal_delay=[zeros(1,n_delay),Phase_code_signal(1:length(Phase_code_signal)-n_delay)];
    fd=2*v_target*fc/c;           %计算多普勒频率
    t=1/fs:1/fs:length(Phase_code_signal_delay)/fs;  
    duppler=exp(1i*2*pi*fd*t);    %产生多普勒频移因子
    Phase_code_signal_return=A_target*Phase_code_signal_delay.*duppler;
end

 function [Phase_code_signal]=Phase_code_modulation(u,fd,fc,fs,Duty_Ratio,n_add)  
      %u为m序列   fd为码频    fc为载频   fs为采样频率
      %产生伪随机相位编码信号
   t = 0:1/fs:(1/fd*length(u)/Duty_Ratio*n_add-1/fs);      % 一个码元所持续的时间内的采样点时刻
   code_signal= exp(1i*2*pi*fc*t);
   N=fs/fd;
   ck=[1, N*length(u)];
   for i=1:length(u)
       for j=1:N
           ck((i-1)*N+j)=u(i);
       end
   end
   ck_zeros=[ck,zeros(1,length(ck)*(1-Duty_Ratio)*10)];
    m_cpi= repmat(ck_zeros,1,n_add);
    Phase_code_signal=  m_cpi.*code_signal;
%    Phase_code_signal= ck_zeros.*code_signal;
%    Phase_code_signal=repmat(Phase_code_signal,1,n_add);

 end
 
 function [signal_demodulation]=demodulation(signal,fc,fs)
  %对回波信号 进行解调
   [~,b] = find(signal~=0);  
   n=length(signal);
   t=[zeros(1,b(1)-1),0:1/fs:(n-b)/fs];
   demodulate= exp(-1i*2*pi*fc*t);
   signal_demodulation=signal.*demodulate;
 end
 
 function [s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR_Preset,xnoise,s_length,N_filter)
    %  add_SNR_xnoise在信号中加入固定信噪比噪声
    %  s_Preset为待加入噪声的信号    SNR_Preset为信噪比  A_xnoise为噪声的标准差
    %  s_length为单周期信号的长度
    %  s为信号     s_add_noise 加入噪声后的信号   SNR_in输入信噪比    xnoise为噪声 
    %  N_filter为滤波器的阶数
    [~,b] = find(s_Preset~=0);    
%     xnoise = A_xnoise*randn(size(s_Preset))+1i*A_xnoise*randn(size(s_Preset));     %产生均值为0，方差为A_xnoise平方的高斯白噪声
    Pn_add=sum(xnoise .*conj(xnoise ))/length(xnoise);                                % 计算加入噪声的平均功率 
    Ps_add=sum(s_Preset(b(1)+N_filter:b(1)+s_length-1).*conj(s_Preset(b(1)+N_filter:b(1)+s_length-1)))/(s_length-N_filter);            % 计算单位信号的平均功率 
    A_1= sqrt((10^(SNR_Preset/10))*(Pn_add/Ps_add));     %计算对应信噪比下信号的幅度
    s=A_1.*s_Preset;                                    %得到对应信噪比下的信号
    s_add_noise=s+xnoise;                            % 加入噪声的信号
    Ps_add_out=sum(s(b(1)+N_filter:b(1)+s_length-1).*conj(s(b(1)+N_filter:b(1)+s_length-1)))/(s_length-N_filter);
    SNR_in=10*log10(Ps_add_out/Pn_add);                   %输入信噪比
 end

  function [s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise_twotarget(s_Preset_1,s_Preset_2,SNR_Preset,A_xnoise,s_length)
    %  add_SNR_xnoise在双回波信号中加入固定信噪比噪声
    %  s_Preset为待加入噪声的信号    SNR_Preset为信噪比  A_xnoise为噪声的标准差
    %  s_length为单周期信号的长度
    %  s为信号     s_add_noise 加入噪声后的信号   SNR_in输入信噪比    xnoise为噪声 
    [~,b] = find(s_Preset_1~=0);    
    xnoise_1 = A_xnoise*randn(size(s_Preset_1))+ 1i*A_xnoise*randn(size(s_Preset_1));     %产生均值为0，方差为A_xnoise平方的高斯白噪声
    xnoise_2 = A_xnoise*randn(size(s_Preset_2))+ 1i*A_xnoise*randn(size(s_Preset_1));     %产生均值为0，方差为A_xnoise平方的高斯白噪声
    if  length(xnoise_1 )>=length(xnoise_2 )
        xnoise=xnoise_1;
    else
        xnoise=xnoise_2;
    end
    Pn_add=sum(xnoise .*conj(xnoise ))/length(xnoise);                                % 计算加入噪声的平均功率 
    Ps_add_1=sum(s_Preset_1(b(1):b(1)+s_length-1).*conj(s_Preset_1(b(1):b(1)+s_length-1)))/(s_length);            % 计算单位信号的平均功率 
    [~,a] = find(s_Preset_2~=0); 
    Ps_add_2=sum(s_Preset_2(a(1):a(1)+s_length-1).*conj(s_Preset_2(a(1):a(1)+s_length-1)))/(s_length);            % 计算单位信号的平均功率 
    Ps_add=(Ps_add_1+Ps_add_2)/2;
    A_1= sqrt((10^(SNR_Preset/10))*(Pn_add/Ps_add));     %计算对应信噪比下信号的幅度
     %实现两个不同长度的信号相加
      r1=length(s_Preset_1);
      r2=length(s_Preset_2);
      s_sum=zeros(1,max(r1,r2));
      s_sum(1,1:r1)=s_Preset_1;
     
      %合成回波信号
      s_sum(1,1:r2)=s_sum(1,1:r2)+s_Preset_2;
      s=A_1*s_sum;
      s_add_noise=s+xnoise;                            % 加入噪声的信号
      Ps_add_out=A_1^2*(Ps_add_1+Ps_add_2)/2;
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
   SK_out=SK.*HK;                %LFM信号与MF做FFT然后相乘
   s_mf_out=ifft(SK_out,N_s);    %得到脉压后的信号
  
   N_n=2^ceil(log2(N3+N2-1));
   NK=fft(xnoise_in,N_n);
   HK_n=fft(h_in,N_n);
   NK_out=NK.*HK_n; 
   xnoise_mf_out_temp=ifft(NK_out,N_n);    %得到脉压后的噪声
   xnoise_mf_out=xnoise_mf_out_temp(N_n/2-N_n/2*(N_s/N_n)+1:N_n/2+N_n/2*(N_s/N_n));  %噪声序列截断
  end

function  [SNR_out,p_s_max,p_n]=SNR_max_out(s,xnoise)
  %求输出信噪比及信号峰值功率、噪声平均功率
  p_s_max=10*log10(max(s)*conj(max(s)));
  p_n=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
  SNR_out =p_s_max- p_n;
end



function  [SNR_out,p_s,p_n]=SNR_caj_time(s,xnoise)
  %求信噪比及功率
  p_s=10*log10(sum(s.*conj(s))/length(s));
  p_n=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
  SNR_out =p_s- p_n;
end

function  [signal_out]=filter_compose(signal,A_filter,b_filter,s_length,Duty_Ratio,N_add)
    [~,b] = find(signal~=0); 
      signal_temp=zeros(N_add,s_length);
     for i=1:N_add
      signal_temp(i,:)=[filter(b_filter,A_filter,signal(b(1)+(i-1)*s_length:(i-1)*s_length+b(1)+s_length*Duty_Ratio-1)),zeros(1,s_length*(1-Duty_Ratio))];
     end
     
      signal_out=zeros(1,b(1)-1);
     for i=1:N_add
         signal_out=[signal_out,signal_temp(i,:)];
     end
 
%      signal_out=[zeros(1,b-1), signal_out];
     signal_out=signal_out(1:length(signal));
end

function [s,s_add_noise,SNR_in,xnoise]=add_SNR(s_Preset,SNR_Preset,xnoise)
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


