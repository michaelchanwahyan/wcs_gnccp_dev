function Y = search_min_vertex(G, L)
    [M, N] = size(G);
    G_ = G;
    trialNum = N; %M*N;
    candi_minVal = zeros(1,trialNum);
    Y_candi = zeros(M*trialNum, N);
    for trialCnt = 1 : trialNum
        Y = zeros(M, N);
        G = G_;
        for mvr = 1 : trialCnt-1 % minimal val removal
            [~,minIdxPos] = min(G(:));
            m_ = mod(minIdxPos-1,M)+1;
            n_ = floor((minIdxPos-1)/M) + 1;
            G(m_,n_) = 999;
        end
        for l_ = 1 : L
            %[~,idx] = sort(G(:));
            %[~,minIdxPos] = min(idx);
            minIdxPos = find(G == min(G(:)));
            m_ = mod(minIdxPos-1,M)+1;
            n_ = floor((minIdxPos-1)/M) + 1;
            Y(m_,n_) = 1;
            G(m_,:) = 999;
            G(:,n_) = 999;
        end
        candi_minVal(trialCnt) = G_(:)' * Y(:);
        Y_candi( (trialCnt-1)*M+1 : trialCnt*M , :) = Y;
        %disp(candi_minVal);
    end
    [~,min_idx] = min(candi_minVal);
    Y = Y_candi( (min_idx-1)*M+1 : min_idx*M , : );
end