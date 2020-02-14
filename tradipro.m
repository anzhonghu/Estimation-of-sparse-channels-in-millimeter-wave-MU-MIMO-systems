H1_e =  U1' * (H1 + (1/sqrt(rho)) * randn(n_arr(1, 1), K));
H2_e =  U2' * (H2 + (1/sqrt(rho)) * randn(n_arr(2, 1), K));
H3_e =  U3' * (H3 + (1/sqrt(rho)) * randn(n_arr(3, 1), K));
diff = 1e4 * ones(n_arr(1, 1), 1);
D1 = zeros(Mb, 1);
D2 = zeros(Mb, 1);
D3 = zeros(Mb, 1);
[~, uhe1in] = sort(abs(H1_e), 'descend');
for k = 1 : K
    counttemp = 0;
    for n = 1 : n_arr(1, 1)
        if diff(uhe1in(n,k), 1) > 0
            diff(uhe1in(n,k), 1) = 0;
            counttemp = counttemp + 1;
            D1((k-1)*beam_numberforMS+counttemp, 1) = uhe1in(n,k);
        else
        end
        if counttemp < beam_numberforMS
        else
            break;
        end
    end
end
H1eb = H1_e(D1, :);
U1b = U1(:, D1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff = 1e4 * ones(n_arr(2, 1), 1);
[~, uhe2in] = sort(abs(H2_e), 'descend');
for k = 1 : K
    counttemp = 0;
    for n = 1 : n_arr(2, 1)
        if diff(uhe2in(n,k), 1) > 0
            diff(uhe2in(n,k), 1) = 0;
            counttemp = counttemp + 1;
            D2((k-1)*beam_numberforMS+counttemp, 1) = uhe2in(n,k);
        else
        end
        if counttemp < beam_numberforMS
        else
            break;
        end
    end
end
H2eb = H2_e(D2, :);
U2b = U2(:, D2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff = 1e4 * ones(n_arr(3, 1), 1);
[~, uhe3in] = sort(abs(H3_e), 'descend');
for k = 1 : K
    counttemp = 0;
    for n = 1 : n_arr(3, 1)
        if diff(uhe3in(n,k), 1) > 0
            diff(uhe3in(n,k), 1) = 0;
            counttemp = counttemp + 1;
            D3((k-1)*beam_numberforMS+counttemp, 1) = uhe3in(n,k);
        else
        end
        if counttemp < beam_numberforMS
        else
            break;
        end
    end
end
H3eb = H3_e(D3, :);
U3b = U3(:, D3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = H1eb / (H1eb' * H1eb + (1/rho)*eye(K));
GI = G' * U1b' * H1;
for k = 1 : K
    interf = 0;
    for kk = 1 : K
        if kk == k
        else
            interf = interf + abs(GI(k, kk))^2;
        end
    end
    interf = rho * interf + (norm(GI(:,k)))^2;
    sinr = rho * abs(GI(k, k))^2  / interf;
    capacity(signalpo_n, 4) = capacity(signalpo_n, 4) + log2(1 + sinr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = H2eb / (H2eb' * H2eb + (1/rho)*eye(K));
GI = G' * U2b' * H2;
for k = 1 : K
    interf = 0;
    for kk = 1 : K
        if kk == k
        else
            interf = interf + abs(GI(k, kk))^2;
        end
    end
    interf = rho * interf + (norm(GI(:,k)))^2;
    sinr = rho * abs(GI(k, k))^2  / interf;
    capacity(signalpo_n, 5) = capacity(signalpo_n, 5) + log2(1 + sinr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = H3eb / (H3eb' * H3eb + (1/rho)*eye(K));
GI = G' * U3b' * H3;
for k = 1 : K
    interf = 0;
    for kk = 1 : K
        if kk == k
        else
            interf = interf + abs(GI(k, kk))^2;
        end
    end
    interf = rho * interf + (norm(GI(:,k)))^2;
    sinr = rho * abs(GI(k, k))^2  / interf;
    capacity(signalpo_n, 6) = capacity(signalpo_n, 6) + log2(1 + sinr);
end