%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  分析讨论信号频点位置和FFT后输出信噪比的关系，说明频点位置造成的最大信噪比损失是多少？
%   对于那些频率分辨率整数倍位置上频点通过加不同的窗分别有多少的信噪比增益损失
%  2021/12/1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  mean(gain_OUT_s)
%%%%%%%%%%%%%%%%% 在两个频点之间分割成若干个频点依次求FFT 的增益 %%%%%%%%%%%%%%%5
times=1000;      %$假定运行100次，按照你的要求修改
Gain_1_s= zeros(1,times) ;
% Gain_2_s= zeros(1,times) ;
D=linspace(0,1,times);
for i_avr=1:times    
  FFT_Gain_noise_B_boxcar  %你要运行的代码
  Gain_1_s(i_avr)=Gain_FFT;
%   Gain_2_s(i_avr)=Gain_2;
end

 f=linspace(1.25e6,1.25e6+fs/N,length(D));
 plot(f,abs(Gain_1_s),'linewidth',1);

%%%%%%%%%%%%%%%%5 求两个分辨率点中点的增益（此时 损耗最大）  %%%%%%%%%%%%%%%%5
% i=1;
% D=0.5;
% FFT_Gain_noise_B_boxcar 
% Gain_FFT_min=Gain_FFT;
% 
% max(Gain_1_s)-min(Gain_1_s)；

