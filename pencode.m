function y=pencode(u,FZlookup,crc_size,bitreversedindices,F_kron_n) 





x = FZlookup;
switch crc_size
        case 0
            crc_code = [];
        case 8
            L=length(u);
            crc_gen=[1 0 0 0 0 0 1 1 1] ;        % CRC generator sequence
            left_shift=[1 0 0 0 0 0 0 0 0];
            a=conv(u,left_shift);        % ×óÒÆ24Î»£¬¼´³ËÒÔ2^24
            for i=1:L                               % Ä£2³ý
                if a(i)==1
                    a(i:i+8)=xor(a(i:i+8),crc_gen);
                end
            end
            crc_code=a(L+1:L+8);                   % CRCÂë


        case 16
            L=length(u);   
            crc_gen=[1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1];        % CRC generator sequence
            left_shift=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            a=conv(u,left_shift);               % ×óÒÆ16Î»£¬¼´³ËÒÔ2^16
            for i=1:L                               % Ä£2³ý
                if a(i)==1
                    a(i:i+16)=xor(a(i:i+16),crc_gen);
                end
            end
            crc_code=a(L+1:L+16);                   % CRCÂë

        case 32
            L=length(u);
            crc_gen=[1 0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 1 0 0 0 1 1 1 0 1 1 0 1 1 0 1 1 1];        % CRC generator sequence
            left_shift=[1 zeros(1,32)];
            a=conv(u,left_shift);        % ×óÒÆ16Î»£¬¼´³ËÒÔ2^16
            for i=1:L                               % Ä£2³ý
                if a(i)==1
                    a(i:i+32)=xor(a(i:i+32),crc_gen);
                end
            end
            crc_code=a(L+1:L+32);   % CRCÂë
end
    
u = [u crc_code];

x (x == -1) = u; % -1's will get replaced by message bits below
x = x(bitreversedindices+1);
y = mod(x*F_kron_n,2);

end