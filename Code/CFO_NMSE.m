function MSE=CFO_NMSE(SNR,Mode,range)
N = 16;

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
        end
    end
end
IDFnT0=DFnT0';

if Mode==0
    Block_Num=1000; %Block Number
    M=4; %Constellation
    w0 = randi([-range,range])/10000;
    L = 2;
    h=(1/sqrt(2*(L+1)))*(randn(1,L+1)+1i*randn(1,L+1));
    N=16;
    K=12;
    Tzp = zeros(N,K);
    Tzp(1:K,1:K) = eye(K);
    Rss=eye(K);    
    %% Matrices
    IFFT=zeros(N);
    for a=1:N
        for b=1:N
            IFFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N);
        end 
    end
    IFFT=IFFT*1/sqrt(N);
    FFT=conj(IFFT);
    Df=zeros(N); 
    for n=1:N
        Df(n,n)=exp(-1i*pi*w0*(n-1));
    end
    
    D=diag(fft(h,16));
    Dk = D(1:K,1:K);
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
    for a=1:Block_Num
        Symbols2(:,:,a)=Df*IDFnT0*H*Tzp*Symbols(:,:,a)+Noise(:,:,a);
    end
    %% CFO Synchronization

    
    Ryy=zeros(N);
    for c=1:Block_Num
        Ryy=Ryy+Symbols2(:,:,c)*Symbols2(:,:,c)';
    end
    Ryy=Ryy/Block_Num;
    J=zeros(20001,1);
    Index=0;
    for w=-1:0.0001:1
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
    
    [M,I] = min(J);

    w_est = (I-10001)/10000;
    
    %% NMSE 
    MSE = (w0-w_est)^2;
else
    N=16;   %Number of Subcarrier
    L=3;    %Channel Length
    Block_Num=1000; %Block Number
    M=4;   %Modulation QAM
    C=4;    %Len Cyclic Prefix, for the CFO estimation scheme, C>L
    W1=randi([-range,range])/10000;  %Subcarrier Frequency Offset 
    W = W1*N/2;
    P=N+C;
    %% Matrix Initialization
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
end