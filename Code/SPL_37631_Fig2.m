loop_Num=1e4;

step=4;
Start=10;
End=30;
total=zeros(1,length(Start:step:End),4);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,length(Start:step:End),4);

for SNRdB=Start:step:End
    disp(SNRdB);
    SNR=10^(SNRdB/10);
    for loop=1:loop_Num
        for Mode=0:3
            [BER]=Fig_2_Final_2(SNR,Mode);
            ratio(1,SNRdB/step-Start/step+1,Mode+1)=BER;
            total(1,SNRdB/step-Start/step+1,Mode+1)=total(1,SNRdB/step-Start/step+1,Mode+1)+ratio(1,SNRdB/step-Start/step+1,Mode+1);
        end
    end
end

total=total/loop_Num;

total2=zeros(1,length(Start:step:End),2);  %preallocating for Speed, SNR from 0 to 20
ratio2=zeros(1,length(Start:step:End),2);

for SNRdB=Start:step:End
    disp(SNRdB);
    SNR=10^(SNRdB/10);
    for loop=1:loop_Num
        for Mode=0
            [BER]=Fig_2_Final_1(SNR,Mode);
            ratio2(1,SNRdB/step-Start/step+1,Mode+1)=BER;
            total2(1,SNRdB/step-Start/step+1,Mode+1)=total2(1,SNRdB/step-Start/step+1,Mode+1)+ratio2(1,SNRdB/step-Start/step+1,Mode+1);
        end
    end
end

total2=total2/loop_Num;

hold on;box on;grid on;
plot(Start:step:End,total(:,:,2),'r--x');
plot(Start:step:End,total(:,:,3),'r:*','linewidth',1.5);
plot(Start:step:End,total(:,:,4),'r-o');
% plot(Start:step:End,total2(:,:,2),'b:*','linewidth',1.5);
plot(Start:step:End,total2(:,:,1),'b-o');
plot(Start:step:End,total(:,:,1),'k-o');

legend('OCDM-NSC:ZF','OCDM-NSC:MMSE','OCDM-NSC:MLE','OCDM[5]:MLE','OFDM[10]:MLE')

set(gca,'Yscale','log');
ylim([1e-8 1e-1]);
xlabel('SNR (dB)');
ylabel('BER');