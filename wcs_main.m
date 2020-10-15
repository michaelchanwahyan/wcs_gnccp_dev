close all
clear
clc

rng(4,'philox');

% currdir = pwd;
% cd /Users/pikachu/Documents/MATLAB/cvx
% cvx_startup
% cd(currdir);

%%
M = 9;
N = 10;

%L = ceil(min(M,N) * 0.5);
L = min(M,N);

VH = 10 * (rand(N,2) - 0.5);

r_ = randintrlv(1:N,1);
r_ = r_(1:M);
r_ = sort(r_);
X_true = zeros(M,N);
for i = 1 : M
    X_true(i,r_(i)) = 1;
end
VG = VH( r_ , : );


%%
AG = get_affinity( VG , 5 );
AH = get_affinity( VH , 5 );

figure;
    subplot(2,2,[1,2]); hold on;
        plot(VH(:,1), VH(:,2), 'bo', 'MarkerSize', 10);
        plot(VG(:,1), VG(:,2), 'rx', 'MarkerSize', 10);
        legend('H', 'G');
        xlim([-5,5]);
        ylim([-5,5]);
        axis equal;
        grid on;
    subplot(2,2,3);
        %pcolor(AG(M:-1:1,:));
        imshow(kron(AG/max(AG(:)), ones(100)));
        title('AG');
        axis equal;
    subplot(2,2,4);
        %pcolor(AH(N:-1:1,:)); axis equal;
        imshow(kron(AH/max(AH(:)), ones(100)));
        title('AH');
        axis equal;

%%
t_intvl = 0.001;
t_arr = 1 : -t_intvl : -1 + t_intvl;
iterNum = length(t_arr);
obj_GNCCP = zeros(1,iterNum);
obj_wcs = zeros(1,iterNum);

% decision variables
X_curr = ones(M,N) * (L/(M*N));

% constants
ONES_NN = ones(N,N);
ONES_MN = ones(M,N);
AG2 = AG.*AG;
AG2T_PLUS_AG2 = AG2' + AG2;
for itr = 1 : iterNum
    t = datetime('now');
    fprintf('[ %s ] @ iter = %d / %d\n', t,  itr, iterNum);
    t = t_arr(itr);
    alpha_t = 1 - abs(t);
    
    % Frank Wolfe Update
    iterNum_fw = 1;
    obj_fw = zeros(1,iterNum_fw);
    
    % search vertex points
        for itr_fw = 1 : iterNum_fw
            % GRAD_J takes nabla H1(X) implementation
                G1 = AG2T_PLUS_AG2 * X_curr * ONES_NN;
                G2 = -2 * (AG' * X_curr * AH + AG * X_curr * AH');
                XtX = X_curr' * X_curr;
                G3 = 2 * (X_curr * AH * XtX * AH' + X_curr * AH' * XtX * AH);
                G = G1 + G2 + G3;
            %{
            Y = zeros(M,N);
            for l_ = 1 : L
                %[~,argidx] = min(G(:));
                argidx = find(G == min(G(:)));
                m_ = mod(argidx,M-1)+1;
                n_ = floor(argidx/M) + 1;
                Y(m_,n_) = 1;
                G(m_,:) = inf;
                G(:,n_) = inf;
            end
            %}
            Y = search_min_vertex(G,L);
            %{
            cvx_begin quiet
            cvx_solver MOSEK
                variable Y(M,N) binary
                minimize( trace(G'*Y) )
                subject to
                    sum(Y,1) <= 1
                    sum(Y,2) <= 1
                    sum(Y(:)) == L
            cvx_end
            Y = round(full(Y));
            %}
            D = Y - X_curr;
            A_ = (D * ONES_NN * D') .* AG - D * AH * D';
            B_ = (X_curr * ONES_NN * D' + ...
                  D * ONES_NN * X_curr') .* AG - ...
                  X_curr * AH * D' - D * AH * X_curr';
            C_ = (X_curr * ONES_NN * X_curr') .* AG - ...
                  X_curr * AH * X_curr';
            a_ = alpha_t * sum(A_(:).^2);
            b_ = 2 * alpha_t * (A_(:)'*B_(:));
            c_ = alpha_t * ( 2 * A_(:)'*C_(:) + sum(B_(:).^2) ) + t * sum(D(:).^2);
            d_ = 2 * ( alpha_t * (B_(:)'*C_(:)) + t * (X_curr(:)'*D(:)) );
            l_ = 0 : 0.001 : 1;
            lambda_ = l_; f = lambda_ * d_;
            lambda_ = lambda_ .* l_; f = f + lambda_ * c_;
            lambda_ = lambda_ .* l_; f = f + lambda_ * b_;
            lambda_ = lambda_ .* l_; f = f + lambda_ * a_;
            % figure; plot(l_, f);
            [~,argmin] = min(f);
            X_prev = X_curr;
            X_curr = X_curr + l_(argmin) * D ;
            %X_curr = awgn(X_curr,100,'measured') + l_(argmin) * D ;
            obj_fw(itr_fw) = gnccp_obj(X_curr, t, AG, AH);
            if (argmin == 1)
                % disp(obj_fw);
                break;
            end
        end
    
    obj_GNCCP(itr) = gnccp_obj(X_curr, t, AG, AH); % Jt(X)
    obj_wcs(itr) = wcs_obj(X_curr, AG, AH);
    %%{
    if (itr > 1 && obj_wcs(itr) > obj_wcs(itr-1))
        X_curr = X_prev;
        obj_GNCCP = obj_GNCCP(1:itr-1);
        obj_wcs = obj_wcs(1:itr-1);
        break;
    end
    %%}
    % fprintf('    obj val = %f\n', obj_GNCCP(itr));
end

X = search_max_vertex(X_curr,L);
%{
cvx_begin quiet
cvx_solver MOSEK
    variable X(M,N) binary
    maximize( trace(X_curr'*X) )
    subject to
        sum(X,1) <= 1
        sum(X,2) <= 1
        sum(X(:)) == L
cvx_end
X = round(full(X));
%}

fprintf('projected soln gnccp_obj = %f\n', gnccp_obj(X, t, AG, AH));
fprintf('projected soln wcs_obj = %f\n', wcs_obj(X, AG, AH));

figure;
plot(obj_GNCCP); hold on;
plot(obj_wcs);
xlabel('it');
title('objective value');
legend('GNCCP', 'WCS');

disp([VG,X_true*VH,X*VH]);

xd = X - X_true;
disp(sum(xd(:) ~= 0) / 2);
figure;
subplot(2,2,1); imshow(X); title('est.');
subplot(2,2,2); imshow(X_true); title('true');
UAG = (X*ONES_NN*X').*AG;
subplot(2,2,3); imshow(kron(UAG/max(UAG(:)),ones(100))); title('U.AG');
XAHXt = X*AH*X';
subplot(2,2,4); imshow(kron(XAHXt/max(XAHXt(:)),ones(100))); title('XA_HX^T');