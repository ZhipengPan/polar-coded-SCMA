function mhat_llr = JIDD(y,polar_N,polar_K,FZlookup,K,V,M,N,CB,N0,h,iter_num,isInterleaver,interleaver,alpha)

    F = zeros(K, V);
		s = [K, M];
		for k = 1:V
			IND = find(CB(:,:,k));
			[I, ~] = ind2sub(s, IND);
			F(unique(I),k) = 1;
		end
       
		LLR = zeros(V, log2(M)*N);
        f = zeros(M, M, M, K,N);
        mhat_llr = zeros(V,polar_K);
        u_llr = zeros(polar_N,V);
        
        for jj = 1:N
            for k = 1:K % resourses
                ind = find(F(k,:)==1); % non-zero elements, paths
                for m1 = 1:M
                    for m2 = 1:M
                        for m3 = 1:M
                            f(m1,m2,m3,k,jj) = -(1/N0)*abs(y(k,jj)-(CB(k,m1,ind(1))*h(k,ind(1),jj)+CB(k,m2,ind(2))*h(k,ind(2),jj)+CB(k,m3,ind(3))*h(k,ind(3),jj)  ))^2;
                        end
                    end
                end
            end
        end
       Ap = 1/M*ones(1,V,M,N);
       Ivg = log(1/M*ones(K, V, M,N));
       Igv = zeros(K, V, M,N);
     
       
       %µü´ú¿ªÊ¼
        for out_loop_iter = 1:iter_num
            for k = 1:K
                ind = find(F(k,:)==1);
                for m1 = 1:M
                    sIgv = zeros(M*M,N);
                    for m2 = 1:M
                        for m3 = 1:M
                            sIgv((m2-1)*M+m3, :) = reshape(f(m1,m2,m3,k, :),1,N)+reshape(Ivg(k,ind(2),m2, :),1,N)+reshape(Ivg(k,ind(3),m3, :),1,N);
                        end
                    end
                    Igv(k,ind(1),m1, :) = log_sum_exp(sIgv);
                end

                for m2 = 1:M
                    sIgv = zeros(M*M, N);
                    for m1 = 1:M
                        for m3 = 1:M
                            sIgv((m1-1)*M+m3,:) = reshape(f(m1,m2,m3,k,:),1,N)+reshape(Ivg(k,ind(1),m1,:),1,N)+reshape(Ivg(k,ind(3),m3, :),1,N);
                        end
                    end
                    Igv(k,ind(2),m2,:) = log_sum_exp(sIgv);
                end

                for m3 = 1:M
                    sIgv = zeros(M*M, N);
                    for m1 = 1:M
                        for m2 = 1:M
                            sIgv((m1-1)*M+m2,:) = reshape(f(m1,m2,m3,k,:),1,N)+reshape(Ivg(k,ind(1),m1,:),1,N)+reshape(Ivg(k,ind(2),m2,:),1,N);
                        end
                    end
                    Igv(k,ind(3),m3,:) = log_sum_exp(sIgv);
                end
            end
            
            Q = zeros(M, V, N);
            for k = 1:V
                ind = find(F(:,k)==1);
                for m = 1:M
                    Q(m,k,:) = reshape(Igv(ind(1),k,m,:),1,N)+reshape(Igv(ind(2),k,m,:),1,N);
                        %Q(m,k,:) = reshape(Igv(ind(1),k,m,:),1,N)+reshape(Igv(ind(2),k,m,:),1,N);
                end
            end
            
            for jj = 1:N
                for k = 1:V
                    LLR(k,2*jj-1) = log((exp(Q(1,k,jj))+exp(Q(2,k,jj)))/((exp(Q(3,k,jj))+exp(Q(4,k,jj)))));
                    LLR(k,2*jj   ) = log((exp(Q(1,k,jj))+exp(Q(3,k,jj)))/((exp(Q(2,k,jj))+exp(Q(4,k,jj)))));
                end
            end
            

            %de-interleaver
            if isInterleaver~=0
                LLR_deinterleaver = zeros(V,polar_N);
                for ii = 1:V
                    LLR_deinterleaver(ii,interleaver(ii,:)) = LLR(ii,:);
                end
            else
                LLR_deinterleaver = LLR;
            end
            %polar decoding
            
            c_llr = zeros(polar_N,V);
            for user = 1:V
                [u_llr(:,user),c_llr(:,user)] = polar_SCAN_decode_alpha(LLR_deinterleaver(user,:),1,alpha,FZlookup,polar_N);
                %mhat_llr(user,:) = u_llr(PCparams.FZlookup==-1,user)';
            end
            
            c_llr = c_llr';
            if isInterleaver~=0
                c_llr_interleaver = zeros(size(c_llr));
                for ii = 1:V
                    c_llr_interleaver(ii,:) = c_llr(ii, interleaver(ii,:));
                end
            else
                c_llr_interleaver = c_llr;
            end
            
            polar2scma = zeros(1,V,M,N);
            for iii = 1:V
                polar2scma(1,iii,1,:) =(exp(c_llr_interleaver(iii,1:2:end))./(exp(c_llr_interleaver(iii,1:2:end))+1)).*(exp(c_llr_interleaver(iii,2:2:end))./(exp(c_llr_interleaver(iii,2:2:end))+1));%00
                polar2scma(1,iii,2,:) =(exp(c_llr_interleaver(iii,1:2:end))./(exp(c_llr_interleaver(iii,1:2:end))+1)).*(1./(exp(c_llr_interleaver(iii,2:2:end))+1));%01
                polar2scma(1,iii,3,:) =(1./(exp(c_llr_interleaver(iii,1:2:end))+1)).*(exp(c_llr_interleaver(iii,2:2:end))./(exp(c_llr_interleaver(iii,2:2:end))+1));%10
                polar2scma(1,iii,4,:) =(1./(exp(c_llr_interleaver(iii,1:2:end))+1)).*(1./(exp(c_llr_interleaver(iii,2:2:end))+1));%11
          
            end
            Ap = polar2scma;
%             for k = 1:V
%                 ind = find(F(:,k)==1);
%                 for n = 1:M
%                     Ivg(ind(1),k,n,:) = log(Ap(1,k,n,:))+Igv(ind(2),k,n,:)-log(sum(exp(Igv(ind(2),k,:,:))));
%                     Ivg(ind(2),k,n,:) = log(Ap(1,k,n,:))+Igv(ind(1),k,n,:)-log(sum(exp(Igv(ind(1),k,:,:))));
%                 end
%             end

            Ivg = log(repmat(Ap,K,1));
            
        end
        
        for user = 1:V
            mhat_llr(user,:) = u_llr(FZlookup==-1,user)';
        end
end