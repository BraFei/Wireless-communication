 clear all; 
%%%%%% Parameter Setting %%%%%%%%%
N_frame=1000;  N_packet=100;
b=2; M=2^b;  % Number of bits per symbol and Modulation order
%  N=input('Please Input the number of antennas: \n');
N=4;
 NT=N;  NR=N; I=eye(NR*2,NR*2);

 
N_pbits = N_frame*NT*b; 
N_tbits = N_pbits*N_packet;
fprintf('====================================================\n');
fprintf(' Pre-MMSE transmission');
fprintf('\n  %d x %d MIMO\n  %d QAM', NT,NR,M);
fprintf('\n  Simulation bits : %d',N_tbits);
fprintf('\n====================================================\n');
SNRdBs = [0:2:20];

for i_SNR=1:length(SNRdBs)
    SNRdB = SNRdBs(i_SNR); 
    noise_var = NT*0.5*10^(-SNRdB/10); 
    sigma = sqrt(noise_var);
    rand('seed',1); randn('seed',1);  
    N_ebits = 0;   
    %%%%%%%%%%%%% Transmitter %%%%%%%%%%%%%%%%%%
    for i_packet=1:N_packet
        X_RTS_hat=zeros(NT,N_frame);
        mod_obj=modem.qammod('M',M,'SymbolOrder','Gray','InputType','bit');
        demod_obj = modem.qamdemod(mod_obj);
        msg_bit = randi([0,1],N_pbits,1); % bit generation
        symbol = modulate(mod_obj,msg_bit).';
%         Scale = modnorm(symbol,'avpow',1); % normalization
        Scale=sq;
        %% 模型转化复数转化为实数模型
        Symbol_nomalized1 = reshape(Scale*symbol,NT,N_frame); 
        Symbol_nomalized=[real(Symbol_nomalized1);imag(Symbol_nomalized1)];
        H1 = randn(NR,NT)+1j*randn(NR,NT);
        H=[real(H1) -imag(H1);imag(H1) real(H1)];
        %%  noise n
        n1 = sigma*complex(randn(NR,1),randn(NR,1));
        n=[real(n1);imag(n1)];
        %% recevied signal
        y=H*Symbol_nomalized +repmat(n,1,N_frame);
        %%
        temp_W = H'*inv(H*H'+noise_var*I); 
        beta = sqrt(2*NT/trace(temp_W*temp_W')); % Eq.(12.17)
        W = beta*temp_W;  
        Tx_signal = W*Symbol_nomalized;  
        %%%%%%%%%%%%% Channel and Noise %%%%%%%%%%%%%        
        Rx_signal = H*Tx_signal+sigma*(randn(NR*2,N_frame));
        %%%%%%%%%%%%%% Receiver %%%%%%%%%%%%%%%%%%%%%
        s_hat = Rx_signal/beta; % Eq.(12.18)
        %% 将实数模型转化为复数方便接下来解码
        s_hat_new=s_hat(1:NT,:)+sqrt(-1)*s_hat(NT+1:end,:);
        Symbol_hat = reshape(s_hat_new/Scale,NT*N_frame,1);
        %% RTS init
        msg_hat = demodulate(demod_obj,Symbol_hat);
        N_ebits = N_ebits + sum(msg_hat~=msg_bit);
    end
        BER(i_SNR) = N_ebits/N_tbits
end
semilogy(SNRdBs,BER,'-k','LineWidth',1); hold on; grid on;
xlabel('SNR[dB]'), ylabel('BER');
legend('Pre-RTS transmission')