%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% 分析多普勒敏感现象
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc 
% clear 
close all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  m_length=127;                     %码长为127
fm=72e6;                          %码频为72MHz
fc=10e9;                          %信号载频为 10GHz
fs=72e6*4;                         %采样率为200MHz
Ts=1/fs;
Duty_Ratio= 0.10;                  %占空比为10%
PRT=m_length*(1/fm)/Duty_Ratio;    %脉冲重复周期
PRF=1/PRT;                        %脉冲重复频率
n_add=round(100e-3/(m_length/fm)*Duty_Ratio);        %相参累计个数
n_add=128;
CPI=n_add*PRT;                       %相参处理时间
SNR_Preset=-10;                       %输入信噪比
A_xnoise =1;                            %噪声标准差
N_PRT=m_length*(fs/fm)/Duty_Ratio;    %一个周期内的采样点数

T_width=m_length/fm;              %信号的时宽
B=fm;
%%%%%%%%%%%%%%%%%%%%%%  产生m序列   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 高斯随机信号 
u = idinput(m_length, 'prbs',[0,1],[-1,1])';   %产生m序列
%%%%%%%%%%%%%%%%%%%%%%%% 产生伪随机相位编码调制信号 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Phase_code_signal=Phase_code_modulation( u,fm,fc,fs,Duty_Ratio,n_add) ;      %产生发射信号
%    plot( real(Phase_code_signal)),   
s_length=length(Phase_code_signal);                        %发射信号的长度

%%%%%%%%%%%%%%%%%%%%%%%%%%%  产生回波信号 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_target=[0,14200];               %目标的距离    R_max=(3e8/2)*PRT*(1-Duty_Ratio)=14288  R_min=(3e8/2)*PRT*Duty_Ratio=1587
v_target=[30,100];                    %目标的速度    c_min=3e8*PRF/2/fc=141.7
A_target=[1,0.3];                  %目标回波信号的幅值衰减
fd=linspace(0,1/(127/fm)*5,300);
mf_max=zeros(1,length(fd));          %脉压后的主瓣峰值
mf_max_2=zeros(1,length(fd));        %脉压后的主瓣峰值
mf_sec=zeros(1,length(fd));          %脉压后的旁瓣峰值
max_sec=zeros(1,length(fd));         %主旁瓣比值
b_max_zhuban=zeros(1,length(fd));  %储存脉压后最大值的位置。
for i=1:length(fd)

   Phase_code_signal_return_1=Phase_code_signal_return_modulation(Phase_code_signal,fs,fc,Duty_Ratio,R_target(1),v_target(1),A_target(1),fd(i));
   h=conj(fliplr(Phase_code_signal(1:m_length*round(fs/fm))));   %产生匹配滤波器而对冲激响应 
   mf_out=conv(Phase_code_signal_return_1(1:m_length*round(fs/fm)),h);
%    [SLL,S_sec]=SLL_CAL(mf_out);
%   mf_sec(i)=S_sec;
   [mf_max(i),b_max_zhuban(i)]=max(abs(mf_out));
   mf_max_2(i)=abs(mf_out(b_max_zhuban(1)));
%    max_sec(i)=mf_max(i)/mf_sec(i);
end

%  plot(fd,20*log10(abs(mf_max)),'linewidth',2),title('多普勒敏感现象')
%  plot(fd,20*log10(abs(mf_max_2)),'linewidth',2),title('多普勒敏感现象')
%  %保持主瓣位置相同下的多普勒容限。



function  [Phase_code_signal_return]=Phase_code_signal_return_modulation(Phase_code_signal,fs,fc,Duty_Ratio,R_target,v_target,A_target,fd)
    %产生雷达回波信号
    %对发射信号延时，加多普勒频率
    c=3e8;    %光速
    n_delay=round(2*R_target*fs/c);    %由于距离引起的延时
    Phase_code_signal_delay=[zeros(1,n_delay),Phase_code_signal(1:length(Phase_code_signal)-n_delay)];
%     fd=2*v_target*fc/c;           %计算多普勒频率
    t=1/fs:1/fs:length(Phase_code_signal_delay)/fs;  
    duppler=exp(1i*2*pi*fd*t);    %产生多普勒频移因子
    Phase_code_signal_return=A_target*Phase_code_signal_delay.*duppler;
end

 function [Phase_code_signal]=Phase_code_modulation(u,fd,fc,fs,Duty_Ratio,n_add)     %产生一个周期内的伪随机相位编码信号
   %u为m序列   fd为码频    fc为载频   fs为采样频率
   
   t = 0:1/fs:(1/fd*length(u)/Duty_Ratio-1/fs);      % 一个码元所持续的时间内的采样点时刻
   code_signal= exp(1i*2*pi*fc*t);
   N=fs/fd;
   ck=[1, N*length(u)];
   for i=1:length(u)
       for j=1:N
           ck((i-1)*N+j)=u(i);
       end
   end
   ck_zeros=[ck,zeros(1,length(ck)*(1-Duty_Ratio)*10)];
%    m_cpi= repmat(ck_zeos,1,n_add);
   Phase_code_signal= ck_zeros;%.*code_signal;  %不加载频
   Phase_code_signal=repmat(Phase_code_signal,1,n_add);

 end
 
 function [SLL,S_sec]=SLL_CAL(signal)
    %求信号的旁瓣电平
   
%     signal=signal(round(11/30*length(signal)):round(19/30*length(signal)));
    Lmax=diff(sign(diff(abs(signal))))==-2;    %得到极值点
    Lmax=[false,Lmax,false];
    vmax=signal(Lmax);
    [S_max,S_max_index]=max(vmax);             %求出主瓣电平和位置
    S_sec_sum=[vmax(S_max_index-3:S_max_index-1),vmax(S_max_index+1:S_max_index+3)];  %得到主瓣周围的旁瓣
    S_sec=max(S_sec_sum);               %得到最大旁瓣
%     Signal_sort=sort(abs(vmax));    %对序列排序
%     S_max= Signal_sort(length(Signal_sort));
%     S_sec=Signal_sort(length(Signal_sort)-1);
    SLL=20*log10(abs(S_sec)/abs(S_max));
 end