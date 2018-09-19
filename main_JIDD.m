%Joint polar code and SCMA
% Copyright (C) 2017  Zhipeng Pan
%School of Electronic Science£¬ National University of Defense Technology 
%zhipengpan10@163.com
% This program is free software: you can redistribute it and/or modify it 

% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

%the result is save in the Joint_CB?_N_K_sigma_alpha_Gaussian_iter_M_interleaverSize_T_result.txt
%Read codebook of SCMA
tic
EbN0 = 1:5;
polar_N = 256;
polar_K = 128;
polar_n = log2(polar_N);
construction_method = 0;%0-BA,1-MC,2-GA
design_snr_dB = 0;%BA and MC construction method will use it
sigma = 0.9;%GA construction method will use it
crc_size = 0;
[FZlookup,bitreversedindices,F_kron_n] = initPC(polar_N,polar_K,polar_n,construction_method,design_snr_dB,sigma,crc_size); 

alpha = 0.6;
iter_num = 5;
isInterleaver = 1;

load('codebook_6users_4chips_qpsk.mat','CB');

K = size(CB, 1); % number of orthogonal resources
M = size(CB, 2); % number of codewords in each codebook
V = size(CB, 3); % number of users (layers)
%polar initial and encoding



SCAN_ITER_NUM = 1;
N = polar_N/log2(M); %Number of scma symbols of each user
SNR  = EbN0 + 10*log10(polar_K/polar_N*log2(M)*V/K);
N0 = 1./10.^(SNR/10); % Noise variance

Nerr = zeros(1,length(EbN0));
Nbits = zeros(1,length(EbN0));
BER   = zeros(1, length(EbN0));

%maxNumErrs = 10000;
maxNumBits = 1e7; %total numer of bits
minNumBits = 50000;
minNumErrs = 50;




for iter_ebn0 = 1:length(EbN0)

    while ((min(Nerr(:,iter_ebn0)) < minNumErrs) && (Nbits(1,iter_ebn0) < maxNumBits) || (Nbits(1,iter_ebn0) <minNumBits) )%100 010 000
        infobits = randi([0 1],V,polar_K);
        c = zeros(V,polar_N);
        for user = 1:V
            c(user,:) = pencode(infobits(user,:),FZlookup,crc_size,bitreversedindices,F_kron_n); 
        end
        
        if isInterleaver ~= 0
            interleaver = zeros(V,polar_N);
            interleavered_bits = zeros(size(c));
            for ii = 1:V
                interleaver(ii,:) = randperm(polar_N);
                interleavered_bits(ii,:) = c(ii, interleaver(ii,:));
            end
            
        else
            interleavered_bits = c;
        end
        
        temp1 = reshape(interleavered_bits',polar_N*V,1);
        temp2 = reshape(temp1,log2(M),N*V);
        x_temp = bi2de(temp2',log2(M),'left-msb');
        x = reshape(x_temp,N,V);
        x = x';
        %h = 1/sqrt(2)*(randn(K, V, N)+1j*randn(K, V, N)); % Rayleigh channel
        h = ones(K, V, N);
        %h = 1/sqrt(2)*(repmat(randn(1, V, N), K, 1)+1j*repmat(randn(1, V, N), K, 1));
        s = scmaenc(x, CB, h); 
        y = awgn(s, SNR(iter_ebn0),'measured');
        
        %Factor graph calculation
        
        
        
        mhat_llr = JIDD(y,polar_N,polar_K,FZlookup,K,V,M,N,CB,N0(iter_ebn0),h,iter_num,isInterleaver,interleaver,alpha);
        
        %**********************************************************
        llr = reshape(mhat_llr',1,V*polar_K);
        m_reshape = reshape(infobits', 1, polar_K*V);
        m_hat = llr<0;
        err = sum(m_hat~=m_reshape);
        Nerr(iter_ebn0) = Nerr(iter_ebn0) + err;
        Nbits(iter_ebn0) = Nbits(iter_ebn0) + length(m_reshape);     
    end
    

    BER(iter_ebn0) = Nerr(iter_ebn0)/Nbits(iter_ebn0);	
	
    fprintf('EbN0 is %d, have runned %d bits, found %d errors, BER=%.7f \n',EbN0(iter_ebn0),Nbits(iter_ebn0),Nerr(iter_ebn0),BER(iter_ebn0));
   

    
    
end

toc
figure;
semilogy(EbN0,BER,'linewidth',1);
