function [FZlookup,bitreversedindices,F_kron_n] = initPC(N,K,n,construction_method,design_snr_dB,sigma,crc_size) 

    F = [1 0;1 1];
    BB=1;
    for ii=1:n
        BB = kron(BB,F);
    end
    F_kron_n = BB;
    
    bitreversedindices = zeros(1,N);
    for index = 1 : N
        bitreversedindices(index) = bin2dec(wrev(dec2bin(index-1,n)));
    end
    switch construction_method
        case 0
            constructed_code_file_name = sprintf('constructedCode\\PolarCode_block_length_%d_designSNR_%.2fdB_method_BhattaBound.txt',N,design_snr_dB);
        case 1
            constructed_code_file_name = sprintf('constructedCode\\PolarCode_block_length_%d_designSNR_%.2fdB_method_MC.txt',N,design_snr_dB);
        case 2
            constructed_code_file_name = sprintf('constructedCode\\PolarCode_block_length_%d_sigma_%.2f_method_GA.txt',N,sigma);
        otherwise
            error('The range of construction_method is from 0 to 2!');
    end


    %should first use construct_polar_code(n) to construct the polar code
    indices = load(constructed_code_file_name);
    FZlookup = zeros(1,N);
    if crc_size == 0
        FZlookup(indices(1:K)) = -1;
    else
        FZlookup(indices(1:K+crc_size)) = -1;
    end
    
end
