function [u_llr, c_llr] = polar_SCAN_decode_alpha(y_llr,iter_num,alpha,FZlookup,N)
    %初始化PCparams.L 和 PCparams.B

    
    n = log2(N);


    plus_infinity = 1000;
    L = zeros(N,n+1);%left message
    B = zeros(N,n+1);%right message
    L(:,n+1) = y_llr';%initial L
    B(FZlookup==0,1) = plus_infinity;%initial B
    
    
    %主循环
    for ii = 1:iter_num
        for phi = 0:N-1
            [L,B] = updateLLRMap(n,phi,n,L,B);
            if mod(phi,2)~=0
                [L,B] = updateBitMap(n,phi,n,L,B);   
            end
        end
    end
     
    mean_B = mean(abs(B(:,n+1)));
    mean_L = mean(abs(L(:,n+1)));
    %输出最终的左信息u_llr和右信息c_llr
    u_llr = L(:,1)+B(:,1);
    
    c_llr = B(:,n+1)+alpha*mean_B/mean_L*L(:,n+1);
    
    
end