% Author: Paulo Francisco
% Objective: Channel modeling and received signal
% Syntax:
%       [y]=channelModeling (ToAs, AoDs, AoAs, Nt, Nr, N, B, Ns)
% Inputs:
%       ToAs, AoDs, AoAs: Localization ParametersBS 
%       Nt, Nr: Number of Tx and Rx antennas
%       N, B, Ns:   Numer of subcarriers, BW and Number of symbols
%
% Outputs:
%       y - received signal
%       
function [y,x,H]=channelModeling (ToAs, AoDs, AoAs, Nt, Nr, N, B, Ns)
    h=10;           % Gain
    Ts=1/B;         % Sampling period in us
    sigma=0.001;      % Noise standard deviation
    L=size(ToAs);   % Number of paths

    %% Channel Modeling
    array_Tx=getArray(Nt);
    array_Rx=getArray(Nr);    
    H=zeros(Nr,Nt,N);
    for n=1:N
        for l=1:L
            toa=ToAs(l);
            aoa=AoAs(l,:);
            aod=AoDs(l,:);
            H(:,:,n)=H(:,:,n)+h*exp(-1j*2*pi*toa*(n-1)/(N*Ts))*...
            sqrt(Nr)*getResponse(array_Rx,aoa)*...
            sqrt(Nt)*getResponse(array_Tx,aod)';
        end
    end

    %% Generate the received signal 
    y=zeros(Nr,Ns,N);
    x=zeros(Nt,Ns,N);    
    for k=1:Ns
        for n=1:N
            x(:,k,n)=exp(1j*rand(Nt,1)*2*pi);  % randon signal
            signal=H(:,:,n)*x(:,k,n);            
            %w=(randn(size(sinal))*std(sinal)+1j*randn(Nr,1))/db2mag(SNR-10);
            w=sigma/sqrt(2)*(randn(Nr,1)+1j*randn(Nr,1));
            y(:,k,n)=signal+w; 
        end
    end
end