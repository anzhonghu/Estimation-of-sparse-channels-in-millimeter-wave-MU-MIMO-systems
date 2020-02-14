for k = 1 : K
    for n_ar_in = 1 : 3
        h = zeros(n_arr(n_ar_in, 1), 1);
        switch n_ar_in
            case 1
                for nx = 0 : naz(n_ar_in, 1)-1
                    for ny = 0 : nel(n_ar_in, 1)-1
                        n = ny * naz(n_ar_in, 1) + 1 + nx;
                        h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * spatailf1(k, 1) + (ny-0.5*(nel(n_ar_in, 1)-1)) * spatailf1(k, 2)));
                    end
                end
                H1_e(:, k) = h * spatailf1(k, 3);
            case 2
                for nx = 0 : naz(n_ar_in, 1)-1
                    for ny = 0 : nel(n_ar_in, 1)-1
                        n = ny * naz(n_ar_in, 1) + 1 + nx;
                        h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * spatailf2(k, 1) + (ny-0.5*(nel(n_ar_in, 1)-1)) * spatailf2(k, 2)));
                    end
                end
                H2_e(:, k) = h * spatailf2(k, 3);
            case 3
                for nx = 0 : naz(n_ar_in, 1)-1
                    for ny = 0 : nel(n_ar_in, 1)-1
                        n = ny * naz(n_ar_in, 1) + 1 + nx;
                        h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * spatailf3(k, 1) + (ny-0.5*(nel(n_ar_in, 1)-1)) * spatailf3(k, 2)));
                    end
                end
                H3_e(:, k) = h * spatailf3(k, 3);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H1eb = U1b' * H1_e;
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
    capacity(signalpo_n, 1) = capacity(signalpo_n, 1) + log2(1 + sinr);
end
H2eb = U2b' * H2_e;
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
    capacity(signalpo_n, 2) = capacity(signalpo_n, 2) + log2(1 + sinr);
end
H3eb = U3b' * H3_e;
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
    capacity(signalpo_n, 3) = capacity(signalpo_n, 3) + log2(1 + sinr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%upper bound
%         H1eb = U1b' * H1;
%         G = H1eb / (H1eb' * H1eb + (1/rho)*eye(K));
%         GI = G' * H1eb;
%         for k = 1 : K
%             interf = 0;
%             for kk = 1 : K
%                 if kk == k
%                 else
%                     interf = interf + abs(GI(k, kk))^2;
%                 end
%             end
%             interf = rho * interf + (norm(GI(:,k)))^2;
%             sinr = rho * abs(GI(k, k))^2  / interf;
%             capacity(signalpo_n, 7) = capacity(signalpo_n, 7) + log2(1 + sinr);
%         end
%         H1eb = U1b' * H11;
%         G = H1eb / (H1eb' * H1eb + (K/rho) * eye(K));
%         G = G / sqrt(trace(G * G') / rho);
%         GI = H1eb' * G;
%         for k = 1 : K
%             interf = 0;
%             for kk = 1 : K
%                 if kk == k
%                 else
%                     interf = interf + abs(GI(k, kk))^2;
%                 end
%             end
%             interf = interf + 1;
%             sinr =  abs(GI(k, k))^2  / interf;
%             capacity(signalpo_n, 7) = capacity(signalpo_n, 7) + log2(1 + sinr);
%         end
%         H2eb = U2b' * H2;
%         G = H2eb / (H2eb' * H2eb + (1/rho)*eye(K));
%         GI = G' * H2eb;
%         for k = 1 : K
%             interf = 0;
%             for kk = 1 : K
%                 if kk == k
%                 else
%                     interf = interf + abs(GI(k, kk))^2;
%                 end
%             end
%             interf = rho * interf + (norm(GI(:,k)))^2;
%             sinr = rho * abs(GI(k, k))^2  / interf;
%             capacity(signalpo_n, 8) = capacity(signalpo_n, 8) + log2(1 + sinr);
%         end
%         H2eb = U2b' * H21;
%         G = H2eb / (H2eb' * H2eb + (K/rho) * eye(K));
%         G = G / sqrt(trace(G * G') / rho);
%         GI = H2eb' * G;
%         for k = 1 : K
%             interf = 0;
%             for kk = 1 : K
%                 if kk == k
%                 else
%                     interf = interf + abs(GI(k, kk))^2;
%                 end
%             end
%             interf = interf + 1;
%             sinr =  abs(GI(k, k))^2  / interf;
%             capacity(signalpo_n, 8) = capacity(signalpo_n, 8) + log2(1 + sinr);
%         end
%         H3eb = U3b' * H3;
%         G = H3eb / (H3eb' * H3eb + (1/rho)*eye(K));
%         GI = G' * H3eb;
%         for k = 1 : K
%             interf = 0;
%             for kk = 1 : K
%                 if kk == k
%                 else
%                     interf = interf + abs(GI(k, kk))^2;
%                 end
%             end
%             interf = rho * interf + (norm(GI(:,k)))^2;
%             sinr = rho * abs(GI(k, k))^2  / interf;
%             capacity(signalpo_n, 9) = capacity(signalpo_n, 9) + log2(1 + sinr);
%         end
%         H3eb = U3b' * H31;
%         G = H3eb / (H3eb' * H3eb + (K/rho) * eye(K));
%         G = G / sqrt(trace(G * G') / rho);
%         GI = H3eb' * G;
%         for k = 1 : K
%             interf = 0;
%             for kk = 1 : K
%                 if kk == k
%                 else
%                     interf = interf + abs(GI(k, kk))^2;
%                 end
%             end
%             interf = interf + 1;
%             sinr =  abs(GI(k, k))^2  / interf;
%             capacity(signalpo_n, 9) = capacity(signalpo_n, 9) + log2(1 + sinr);
%         end