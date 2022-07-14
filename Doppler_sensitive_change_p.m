clc
close all
clear
p=[63,127,255];
times=length(p);
sll_sum=zeros(times,300);
for i_sum=1:times   
   m_length=p(i_sum) ;
   Doppler_sensitive
   sll_sum(i_sum,:)= mf_max_2;  %mf_max; %
end

figure(4)
hold on
for j=1:times
    plot(fd,20*log10(abs(sll_sum(j,:))),'linewidth',1),title('多普勒敏感现象')
   
end
legend('码长为：63','码长为：127','码长为：255')
hold off


figure(5)
hold on
plot(fd,20*log10(abs(sll_sum(1,:))), '-*' ),title('多普勒敏感现象')
plot(fd,20*log10(abs(sll_sum(2,:))), '->' ),title('多普勒敏感现象')
plot(fd,20*log10(abs(sll_sum(3,:))), '-+' ),title('多普勒敏感现象')
legend('码长为：63','码长为：127','码长为：255')
hold off