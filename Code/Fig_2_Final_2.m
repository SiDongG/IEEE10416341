function [BER]=Fig_2_Final_2(SNR,Mode)
% clear; clc; close all;
Block_Num=100; %Block Number
M=4; %Constellation
L=2; %Channel Order
N=16; %Block S
P=N+L;
% ize
K=12; %Actual Used Sub-carrier 
w0=randi([-100,100])/100; %Normalized CFO
% Mode=3;
% SNR=1000;
Rss = eye(K);
Spacing = [1,3,6,11];
Tds = zeros(N,K);
Tds(1,:) = Rss(1,:);
Tds(3,:) = Rss(2,:);
Tds(5:6,:) = Rss(3:4,:);
Tds(8:11,:) = Rss(5:8,:);
Tds(13:16,:) = Rss(9:12,:);

%% Fresnel Matrices
DFnT0=zeros(N);
if mod(N,2)==0
    for m=1:N
        for n=1:N
            DFnT0(m,n)=sqrt(1/N)*exp(-1i*pi/4)*exp(1i*pi*((m-1)-(n-1))^2/N);
%             DFnT0(m,n)=sqrt(1/N)*exp(1i*pi*((m-1)-(n-1))^2/N);
        end
    end
end
IDFnT0=DFnT0';
%% Channel Generation
h=(1/sqrt(2*(L+1)))*(randn(1,L+1)+1i*randn(1,L+1));
N=16;
K=12;
Tzp = zeros(N,K);
Tzp(1:K,1:K) = eye(K);

c1=1/(2*N);
c2=1/(2*N);

IFFT=zeros(N);
for a=1:N
    for b=1:N
        IFFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N);
    end 
end
IFFT=IFFT*1/sqrt(N);
FFT=conj(IFFT);

Vc1=zeros(N);
Vc2=zeros(N);
for c=1:N
    Vc1(c,c)=exp(-1i*2*pi*c1*(c-1)^2);
    Vc2(c,c)=exp(-1i*2*pi*c2*(c-1)^2);
end

Df=zeros(N); 
for n=1:N
    Df(n,n)=exp(-1i*pi*w0*(n-1));
end

D=diag(fft(h,16));
Dk = D(1:K,1:K);
A = FFT*Vc1'*IFFT*Vc2';
H = IFFT*D*FFT;
%% Modulation
Bits=randi(0:1,[1,K*Block_Num*log2(M)]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
end
Bits3=qammod(Bits2,M)*sqrt(0.5);
Symbols=reshape(Bits3,K,1,Block_Num);
%% Channel and CFO Matrix
nr=randn(N,1,Block_Num);
ni=randn(N,1,Block_Num);
Noise=(1/sqrt(SNR))*(sqrt(2)/2)*(nr+1i*ni);
Symbols2=zeros(N,1,Block_Num);

if Mode == 0 
    for a=1:Block_Num
        Symbols2(:,:,a)=Df*IFFT*D*Tds*Symbols(:,:,a)+Noise(:,:,a);
    end
else
    for a=1:Block_Num
        Symbols2(:,:,a)=Df*IFFT*D*A*Tzp*Symbols(:,:,a)+Noise(:,:,a);
    end
end

%% Calculate Covariance Matrix 
% Ryy=zeros(N);
% for c=1:Block_Num
%     Ryy=Ryy+Symbols2(:,:,c)*Symbols2(:,:,c)';
% end
% Ryy=Ryy/Block_Num;


Ryy=Df*IFFT*Tds*Dk*Rss*Dk'*Tds'*FFT*Df';

J3=zeros(201,1);
Index=0;
for w=-1:0.01:1
    Dff=zeros(N); 
    for n=1:N
        Dff(n,n)=exp(-1i*pi*w*(n-1));
    end
    Index=Index+1;
    for k=1:N-K
        f=zeros(N,1);
        index = Spacing(k);
        for a=1:N
            f(a)=exp(1i*(a-1)*2*pi*index/N);
        end
        J3(Index)=J3(Index)+f'*inv(Dff)*Ryy*Dff*f;
    end
end

%% CFO Synchronization
J=zeros(201,1);
if Mode==0
    
    Ryy=Df*IFFT*Tds*Dk*Rss*Dk'*Tds'*FFT*Df';
    
    Index=0;
    for w=-1:0.01:1
        Dff=zeros(N); 
        for n=1:N
            Dff(n,n)=exp(-1i*pi*w*(n-1));
        end
        Index=Index+1;
        for k=1:N-K
            f=zeros(N,1);
            index = Spacing(k);
            for a=1:N
                f(a)=exp(1i*(a-1)*2*pi*index/N);
            end
            J(Index)=J(Index)+f'*inv(Dff)*Ryy*Dff*f;
        end
    end
else
    Ryy=Df*IDFnT0*H*Tzp*Rss*Tzp'*H'*DFnT0*Df';

    J=zeros(201,1);
    Index=0;
    for w=-1:0.01:1
        Dff=zeros(N); 
        for n=1:N
            Dff(n,n)=exp(-1i*pi*w*(n-1));
        end
        Index=Index+1;
        for k=K+L+1:N
            lns=DFnT0(k,:);
            J(Index)=J(Index)+abs(lns*inv(Dff)*Ryy*Dff*lns');
        end
    end
end
Index=find(J==min(J));
Est_w=1-0.01*Index+0.01;
Est_w= -Est_w;
%% CFO Compensation
Dff=zeros(N); 
for n=1:N
    Dff(n,n)=exp(1i*pi*Est_w*(n-1));
end
Symbols3=zeros(size(Symbols2));
for count=1:Block_Num
    Symbols3(:,:,count)=Dff*Symbols2(:,:,count);
end
%% Equalization 
Symbols5=zeros(size(Symbols));

% OFDM ML, OCDM-NSC: ZF MMSE ML
if Mode==0
    for count=1:Block_Num
        Symbols5(:,:,count) = qam_sphere_decoder(IFFT*D*Tds,Symbols3(:,:,count),M,Symbols(:,:,count),K);
    end
elseif Mode==1
    for count=1:Block_Num
        Symbols5(:,:,count) = pinv(IFFT*D*A*Tzp)*Symbols3(:,:,count);
    end
elseif Mode==2
    for count=1:Block_Num
        B = IFFT*D*A*Tzp;
        Symbols5(:,:,count) = B'*inv(1/SNR*eye(N)+B*B')*Symbols3(:,:,count);
    end
else
    for count=1:Block_Num
        Symbols5(:,:,count) = qam_sphere_decoder(IFFT*D*A*Tzp,Symbols3(:,:,count),M,Symbols(:,:,count),K);
    end
end

%% Demodulation
if M==4
    Symbols6=qamdemod(Symbols5/sqrt(1/2),M);
end
Bitsre=zeros(1,K*Block_Num*log2(M));
start=1;
for count=1:Block_Num
    for k=1:K
        dec=dec2bin(Symbols6(k,1,count),log2(M));
        for n=1:length(dec)
            Bitsre(start)=str2double(dec(n));
            start=start+1;
        end
    end
end
%% Error Count
BER=sum(Bitsre~=Bits)/length(Bits);