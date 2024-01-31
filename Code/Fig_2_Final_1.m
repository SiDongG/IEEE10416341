function BER=Fig_2_Final_1(SNR,Mode)
% clear; clc; close all;
N=16;   %Number of Subcarrier
L=3;    %Channel Length
Block_Num=10; %Block Number
M=4;   %Modulation QAM
C=8;    %Len Cyclic Prefix, for the CFO estimation scheme, C>L
W1=randi([-100,100])/100;  %Subcarrier Frequency Offset 
W = W1*N/2;
P=N+C;
% W = W*NP/2;
%     SNR = 10000;
%     Mode = 1;
%% Matrix Initialization
% DFnT
% DFT
% CFO
IFFT=zeros(N);
for a=1:N
    for b=1:N
        IFFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N);
    end 
end
IFFT=IFFT*1/sqrt(N);
FFT=conj(IFFT);

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

D=zeros(P,P);
for count=1:P
    D(count,count)=exp(1i*W*2*pi*(count-1)/N);
end
%% Symbols Initialization
Bits=randi(0:1,[1,N*Block_Num*log2(M)]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
end
Bits3=qammod(Bits2,M)*sqrt(0.5);
%% Pilot Symbols Assignment
Symbols=reshape(Bits3,N,1,Block_Num);

%% IDFnT
Symbols1=zeros(N,1,Block_Num);

for count=1:Block_Num
    Symbols1(:,:,count)=IDFnT0*Symbols(:,:,count);
end
%% Cyclic Prefix 
S=eye(N);
T=[S(2*N-P+1:N,:);S];
R=[zeros(N,P-N),eye(N)];
Symbols2=zeros(P,1,Block_Num);
for count=1:Block_Num
    Symbols2(:,:,count)=T*Symbols1(:,:,count);
end
%% Construct Channel Matrix
h=(1/sqrt(2*(L)))*(randn(1,L)+1i*randn(1,L));
H0=zeros(P); %Preallocating for speed, H0 is the P by P matrix have the (i,j)th entry h(i-j)
H1=zeros(P); %Preallocating for speed, H1 is the P by P matrix have the (i,j)th entry h(P+i-j)
a=1;
while a<P+1  %generate the channel matrces
    b=1;
    while b<P+1
        if a-b<0 || a-b>L-1
            H0(a,b)=0;
        else
            H0(a,b)=h(a-b+1);
        end
        if P+a-b<0 || P+a-b>L-1
            H1(a,b)=0;
        else
            H1(a,b)=h(P+a-b+1);
        end
        b=b+1;
    end
    a=a+1;
end
%% Construct Noise Matrix 
nr=randn(P,1,Block_Num);
ni=randn(P,1,Block_Num);
Noise=(sqrt(2)/2)*(nr+1i*ni);
%% AWGN Multipath Channel 
Symbols3=zeros(P,1,Block_Num);
for a=1:Block_Num
    Insertion1=Symbols2(:,:,a);
    if a==1
        Insertion2=zeros(P,1);
    else
        Insertion2=Symbols2(:,:,a-1);
    end
    Symbols3(:,:,a)=H0*Insertion1+H1*Insertion2+(1/sqrt(SNR))*Noise(:,:,a);
end
%% Channel Carrier Frequency Offset 
Symbols4=zeros(P,1,Block_Num);
for count=1:Block_Num
    Symbols4(:,:,count)=exp(1i*W*2*pi*P*count)*D*Symbols3(:,:,count);
end
%% CFO Estimation
Sum=0;
for i=1:Block_Num
    for count=(L+1):C
        Sum=Sum+conj(Symbols4(count,1,i))*Symbols4(count+N,1,i);
    end
end
W0=(1/(2*pi))*imag(log(Sum));
MSE=(W0-W)^2;
Dh = diag(fft(h,16));
H = IFFT*Dh*FFT;
%% Simulation
Df=zeros(N); 
for n=1:N
    Df(n,n)=exp(-1i*pi*W*(n-1)/N);
end
Dff=zeros(length(Df)); 
for n=1:N
    Dff(n,n)=exp(1i*pi*W0*(n-1)/N);
end
Symbolss = zeros(size(Symbols));
for count=1:Block_Num
    Symbolss(:,:,count)=Dff*Df*H*IDFnT0*Symbols(:,:,count);
end
%% 
Symbolsr=zeros(size(Symbols));
for count=1:Block_Num
    Symbolsr(:,:,count) = qam_sphere_decoder(H*IDFnT0,Symbolss(:,:,count),M,Symbols(:,:,count),N);
end

%% Demod
if M==4
    Symbolsrr=qamdemod(Symbolsr/sqrt(1/2),M);
end
Bitsre=zeros(1,N*Block_Num*log2(M));
start=1;
for count=1:Block_Num
    for k=1:N
        dec=dec2bin(Symbolsrr(k,1,count),log2(M));
        for n=1:length(dec)
            Bitsre(start)=str2double(dec(n));
            start=start+1;
        end
    end
end

BER=sum(Bitsre~=Bits)/length(Bits);