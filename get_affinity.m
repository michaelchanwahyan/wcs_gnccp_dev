function W = get_affinity(P, theta)

N = size(P,1);
W = zeros(N,N);
%sigm = theta^2;
for i = 1 : N
    p_diff = P(i,:) - P(i:N,:);
    
    %d_ = sum( p_diff.^2 , 2 ) / sigm;
    %W(i,i:N) = exp(-d_);
    
    d_ = sqrt(sum(p_diff.^2, 2));
    
    W(i,i:N) = d_;
    W(i:N,i) = W(i,i:N);
    W(i,i) = 0;
end
