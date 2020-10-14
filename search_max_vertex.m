function Y = search_max_vertex(G, L)
    G = awgn(G,50,'measured');
    [M, N] = size(G);
    G_ = G;
    trialNum = M*N;
    candi_maxVal = zeros(1,trialNum);
    Y_candi = zeros(M*trialNum, N);
    for trialCnt = 1 : trialNum
        Y = zeros(M, N);
        G = G_;
        for mvr = 1 : trialCnt-1 % maximal val removal
            [~,maxIdxPos] = max(G(:));
            m_ = mod(maxIdxPos-1,M)+1;
            n_ = floor((maxIdxPos-1)/M) + 1;
            G(m_,n_) = 999;
        end
        for l_ = 1 : L
            %[~,idx] = sort(G(:));
            %[~,minIdxPos] = min(idx);
            maxIdxPos = find(G == max(G(:)));
            m_ = mod(maxIdxPos-1,M)+1;
            n_ = floor((maxIdxPos-1)/M) + 1;
            Y(m_,n_) = 1;
            G(m_,:) = -999;
            G(:,n_) = -999;
        end
        candi_maxVal(trialCnt) = G_(:)' * Y(:);
        Y_candi( (trialCnt-1)*M+1 : trialCnt*M , :) = Y;
        %disp(candi_minVal);
    end
    [~,max_idx] = max(candi_maxVal);
    Y = Y_candi( (max_idx-1)*M+1 : max_idx*M , : );
end