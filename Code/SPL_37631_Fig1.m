clear; clc; close all;
loop_Num=1000;

step=5;
Start=0;
End=20;
MSE_set=zeros(1,length(Start:step:End),6);
total=zeros(1,length(Start:step:End),6); 
Range=[500,1000,10000];


for SNRdB=Start:step:End
    disp(SNRdB);
    SNR=10^(SNRdB/10);
    for loop=1:loop_Num
        Count=1;
        for Mode=0:1
            for r=1:3
                range=Range(r);
                MSE=CFO_NMSE(SNR,Mode,range);            
                MSE_set(1,SNRdB/step-Start/step+1,Count)=MSE_set(1,SNRdB/step-Start/step+1,Count)+MSE;
                Count=Count+1;
            end
        end
    end
end

MSE_set=MSE_set/loop_Num;
figure()
box on; hold on;grid on

plot(Start:step:End,MSE_set(:,:,6),'k--');
plot(Start:step:End,MSE_set(:,:,1),'r-');
plot(Start:step:End,MSE_set(:,:,1),'ro');
plot(Start:step:End,MSE_set(:,:,2),'r*');
plot(Start:step:End,MSE_set(:,:,3),'r+');

plot(Start:step:End,MSE_set(:,:,1),'r-o');
plot(Start:step:End,MSE_set(:,:,2),'r-*');
plot(Start:step:End,MSE_set(:,:,3),'r-+');
plot(Start:step:End,MSE_set(:,:,4),'k--o');
plot(Start:step:End,MSE_set(:,:,5),'k--*');
plot(Start:step:End,MSE_set(:,:,6),'k--+');

set(gca,'Yscale','log');
ylim([1e-8 1e2]);
xlabel('SNR (dB)');
ylabel('NMSE of CFO');
legend('OCDM [5]','OCDM-NSC','w_0 \in [-0.05\pi,0.05\pi)','w_0 \in [-0.1\pi,0.1\pi)','w_0 \in [-\pi,\pi)')