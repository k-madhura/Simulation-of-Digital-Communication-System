% ENEE 630 Project

clc; clear all; close all;
%%
% Evaluation 1
for k=1:1:100
    [TD1(k), FE1(k) , ber_1(k), fer_1(k)]=func1(0.0625, 0, 100 );   
end
mean_timedelay=mean(TD1,2);
std_timedelay=std(TD1);
%%
%Evaluation 2
for k=1:1:100
    [TD2(k), FE2(k) , ber_2(k), fer_2(k)]=func1(0, 125, 100 );   
end
mean_freqestim=mean(FE2,2);
std_freqestim=std(FE2);
%%
%Evaluation 3
count=1;
for x = -3:0.5:15
    for k=1:1:10
        [TD3(count, k), FE3(count, k), ber3(count, k), fer3(count, k)] = func1(2.5/80, 62.5, x);       
    end
    count=count+1; 
end
ber_avg1=transpose(mean(ber3, 2));
fer_avg1=mean(fer3, 2);
FE_avg=mean(FE3, 2);
TD_avg=mean(TD3, 2);
x1 = -3:0.5:15;
semilogy(x1, ber_avg1);
semilogy(x1, fer_avg1);

%%
function [time_delay, f_0, ber, fer ]= func1(t0, f0, snr1 )

fs= 16000;       % Carrier Frequency(Hz)
Length=0.005;    % 50ms length

% Signal Generation
data= randi([0 1], 1, 664);     % Generation of 664 Random bits
key= zeros(1, 8);               % Generation of key(8 bits)
codeword= ones(1, 128);         % Generation of codeword( 128 bits)
Data1= [codeword key data];     % Signal Frame
N=length(Data1);                % Finding the length of the bit stream

% Modulation
for k=1:1:N
    x(k)=exp(1i*pi*k/4)*(1-2*(Data1(k)));       % Pi/4 BPSK modulation
    %disp(x(k));
end

% Upsampling
x1=upsample(x,16);                              % Upsampling by a factor or 16
%disp(x1);

% Filtering by rrc
rrcFilter = rcosdesign(0.35,6,97);       % Setting filter specifications- alpha=0.35, span=6 ans sps= 97
s = upfirdn(x1,rrcFilter);               % Implementing transmitter FIR filter
s_1=s(292:end-291);                      % Removing the tail of the filter
%disp(s_1);
plot(abs(s))

% Channel
p=zeros(1, 640);                        % Creating a matrix of zeros to pad the signal
s_2=[p s_1 p];                          % Padding the signal with 16*40 zeros on both sides                        
k0=ceil(t0*16);                                % Equivalent delay of k0
s_3=circshift(s_2,k0,2);                 % Shifting the signal by the delay row wise
for k=1:1:length(s_2)
    s(k)=exp((1i*2*pi*f0*k)/(16*fs))*s_3(k);         % Adding frequency uncertainty
end
%disp(s);

R= awgn(s, snr1);                       % Adding average white gaussian noise to the signal
%disp(R);

% Receiver
rrcFilter1 = rcosdesign(0.35,6,97);     % Specifying receiver filter specifications
s1 = upfirdn(R,rrcFilter1);             % Implementing receiver FIR filter
s2=s1(292: end-291);                    % Removing filter tail
plot(abs(s1));

% Downsampling
y1= abs(s2);                % Finding the magnitude of each sample
res= reshape(y1, [16, 880]);            % Reshaping the matrix to 16*880 size to find sum of each 16th sample
res_sum= sum(res, 2);                      % Taking the sum of each row
%disp("sum of 16 cols")
%disp(res_sum);

[M3 I3]=max(res_sum.',[], 2);             % I3 is the offset( has maximum energy)
%     disp(I3);
%     disp(M3);
y_2=downsample(s2, 16, I3-1);             % Downsampling the signal by factor of 16 and offset I3
%disp(y_2);

%DFT
f=[];
for j=1:1:length(y_2)-128
    f(1,j)=max(abs(fft(y_2(1,j:(j+127)), 128, 2)));
end
[M I]=max(f, [], 2);                            % Magnitude and index of the highest frequency

k_0= 41+I;                                                  % Approximating the time delay
time_delay=I/16;
f_0= 125* I;                                                % Frequency estimate
y_7= circshift(y_2,41-I,2);                                  % Adjusting time delay

for k=1:1:length(y_7)
    y(k)=exp(-2*pi*1i*f_0*k/16000).*y_7(k)*exp(-1i*pi*k/4);        % Adjusting the frequency estimate and demodulating
end

y_final=zeros(1,880);

% Representing two symbols with one according to P/4 bpsk
for k=1:1:880
    if (real(y(1,k))<0) && (imag(y(1,k))<0)
        y_final(1,k)=1;
    elseif (real(y(1,k))>0) && (imag(y(1,k))<0)
        y_final(1,k)=0;
    elseif (real(y(1,k))<0) && (imag(y(1,k))>0)
        y_final(1,k)=1;
    else
        y_final(1,k)=0;
    end
end
y_final1=y_final(1,177:end-40);        % Extracting final signal
%disp(y_final1);

%Finding BER and FER
ber1=zeros(1,700);
for l=1:1:664
    if y_final1(1,l)~=data(1,l)
        ber1(1,l)=1;
    else
        ber1(1,l)=0;
    end
end
ber=sum(ber1, 2)/664;
%disp(ber);
fer1=reshape(ber1, [7, 100]);    % Reshaping to find FER with one frame size of 100 bits
fer2=sum(fer1, 2);      % Sum of errors frame wise
fer=sum(fer2, 1)/7;     % Finding FER
%disp(fer);
    
end



