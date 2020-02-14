clear;
naz = [7; 16; 32];
nel = [36; 76; 151];
n_arr = naz .* nel;
K = 10;
Mb = 2 * K;
beam_numberforMS = floor(Mb/K);
Rmin = 10;
Rmax = 100;
rho = 1e10;
rho0 = rho * 1e5;
height = 10;
f = 80*10^9;%1G bandwidth
lambda = 3 * 10^8 / f;
miu = 0.5;
Tc = 500;
N = Tc-K;
pilot_s = [1; 3; 5; K];
Spad = 5;
betam_min = atan(height/Rmax);
betam_max = atan(height/Rmin);
beta_m = 10.3 * pi / 180;%-phi_m in the EL paper
Nite = 1e3;
c_km = zeros(K, beam_numberforMS);
el_in = zeros(Mb, 1);
az_in = zeros(Mb, 1);
capacity = zeros(length(pilot_s),  9);
esti_err = zeros(length(pilot_s),  3);
angle = zeros(1, 2);
U1 = zeros(n_arr(1, 1), n_arr(1, 1));
U2 = zeros(n_arr(2, 1), n_arr(2, 1));
U3 = zeros(n_arr(3, 1), n_arr(3, 1));
Tc_store = [15; 20; 25; 30;  35];
capacity_tc = zeros(length(Tc_store), 2);
for nx = 1 : naz(1, 1)
    for ny = 1 : nel(1, 1)
        angle(1, 2) = (-1+2*ny/nel(1, 1));%el
        angle(1, 1) = (-1+2*nx/naz(1, 1));%az
        n = (ny - 1) * naz(1, 1) + nx;
        for mx = 0 : naz(1, 1)-1
            for my = 0 : nel(1, 1)-1
                m = my * naz(1, 1) + 1 + mx;
                U1(m, n) = exp(-1i * 2 * pi * miu * ((mx-0.5*(naz(1, 1)-1)) * angle(1, 1) + (my-0.5*(nel(1, 1)-1)) * angle(1, 2))) / sqrt(n_arr(1, 1));
            end
        end
    end
end
for nx = 1 : naz(2, 1)
    for ny = 1 : nel(2, 1)
        angle(1, 2) = (-1+2*ny/nel(2, 1));%el
        angle(1, 1) = (-1+2*nx/naz(2, 1));%az
        n = (ny - 1) * naz(2, 1) + nx;
        for mx = 0 : naz(2, 1)-1
            for my = 0 : nel(2, 1)-1
                m = my * naz(2, 1) + 1 + mx;
                U2(m, n) = exp(-1i * 2 * pi * miu * ((mx-0.5*(naz(2, 1)-1)) * angle(1, 1) + (my-0.5*(nel(2, 1)-1)) * angle(1, 2))) / sqrt(n_arr(2, 1));
            end
        end
    end
end
for nx = 1 : naz(3, 1)
    for ny = 1 : nel(3, 1)
        angle(1, 2) = (-1+2*ny/nel(3, 1));%el
        angle(1, 1) = (-1+2*nx/naz(3, 1));%az
        n = (ny - 1) * naz(3, 1) + nx;
        for mx = 0 : naz(3, 1)-1
            for my = 0 : nel(3, 1)-1
                m = my * naz(3, 1) + 1 + mx;
                U3(m, n) = exp(-1i * 2 * pi * miu * ((mx-0.5*(naz(3, 1)-1)) * angle(1, 1) + (my-0.5*(nel(3, 1)-1)) * angle(1, 2))) / sqrt(n_arr(3, 1));
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
in1_l = 42;
in1_u = 132;
in2_l = 220;
in2_u = 648;
in3_l = 900;
in3_u = 2593;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H1 = zeros(n_arr(1, 1), K);
H2 = zeros(n_arr(2, 1), K);
H3 = zeros(n_arr(3, 1), K);
H1_e = zeros(n_arr(1, 1), K);
H2_e = zeros(n_arr(2, 1), K);
H3_e = zeros(n_arr(3, 1), K);
H1eb = zeros(Mb, K);
H2eb = zeros(Mb, K);
H3eb = zeros(Mb, K);
Cor1 = zeros(naz(1, 1), nel(1, 1));
Cor2 = zeros(naz(2, 1), nel(2, 1));
Cor3 = zeros(naz(3, 1), nel(3, 1));
for signalpo_n = 1 : length(pilot_s)
    pilot_n = pilot_s(signalpo_n, 1);
    N_1 = Tc - pilot_n;
    for ii = 1 : Nite
        pos = zeros(K, 3);
        theta = zeros(K, 1);
        phi = zeros(K, 1);
        f_ind1 = zeros(K, 4);
        f_ind2 = zeros(K, 4);
        f_ind3 = zeros(K, 4);
        beta = zeros(K, 2);
        beam_in = zeros(Mb, 1);
        for k = 1 : K
            pos_temp = zeros(1, 2);
            while norm(pos_temp) < Rmin || norm(pos_temp) > Rmax || abs(atan(pos_temp(1, 2) / pos_temp(1, 1))) > pi / 3
                pos_temp(1, 1) = rand(1, 1) * Rmax;
                pos_temp(1, 2) = (rand(1, 1) * 2 - 1) * Rmax;
            end
            pos(k, 1:2) = pos_temp;
            pos(k, 3) = norm(pos_temp);
            pos(k, 3) = norm([pos(k, 3), height]);%distance
            phi(k, 1) = asin(pos_temp(1, 2) / sqrt(pos_temp(1, 2)^2 + (pos_temp(1, 1) * cos(beta_m) + height * sin(beta_m))^2));%az
            theta(k, 1) = asin((pos_temp(1, 1) * sin(beta_m) - height * cos(beta_m)) / pos(k, 3));%el
            for n_ar_in = 1 : 3
                switch n_ar_in
                    case 1
                        f_ind1(k, 1) = (sin(phi(k, 1)) * cos(theta(k, 1)) + 1) * 0.5 * naz(n_ar_in, 1);
                        f_ind1(k, 2) = (sin(theta(k, 1)) + 1) * 0.5 * nel(n_ar_in, 1);
                        f_ind1(k, 3) = (f_ind1(k, 2) - 1) * naz(n_ar_in, 1) + f_ind1(k, 1);
                    case 2
                        f_ind2(k, 1) = (sin(phi(k, 1)) * cos(theta(k, 1)) + 1) * 0.5 * naz(n_ar_in, 1);
                        f_ind2(k, 2) = (sin(theta(k, 1)) + 1) * 0.5 * nel(n_ar_in, 1);
                        f_ind2(k, 3) = (f_ind2(k, 2) - 1) * naz(n_ar_in, 1) + f_ind2(k, 1);
                    case 3
                        f_ind3(k, 1) = (sin(phi(k, 1)) * cos(theta(k, 1)) + 1) * 0.5 * naz(n_ar_in, 1);
                        f_ind3(k, 2) = (sin(theta(k, 1)) + 1) * 0.5 * nel(n_ar_in, 1);
                        f_ind3(k, 3) = (f_ind3(k, 2) - 1) * naz(n_ar_in, 1) + f_ind3(k, 1);
                end
            end
        end
        [~, f_ind1(:, 4)] = sort(f_ind1(:, 3));
        [~, f_ind2(:, 4)] = sort(f_ind2(:, 3));
        [~, f_ind3(:, 4)] = sort(f_ind3(:, 3));
        for kx = 1 : K
            for n_ar_in = 1 : 3
                switch n_ar_in
                    case 1
                        k = f_ind1(kx, 4);
                    case 2
                        k = f_ind2(kx, 4);
                    case 3
                        k = f_ind3(kx, 4);
                end
                Aaz = -min(12 * phi(k, 1)^2 / (70/180*pi)^2, 25);
                Ael = -min(12 * theta(k, 1)^2 / (7/180*pi)^2, 20);
                D0 = -min(-Aaz-Ael, 25);
                D0 = 10^(D0*0.1);
                beta(k, 1) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos(k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi);
                h = zeros(n_arr(n_ar_in, 1), 1);
                for nx = 0 : naz(n_ar_in, 1)-1
                    for ny = 0 : nel(n_ar_in, 1)-1
                        n = ny * naz(n_ar_in, 1) + 1 + nx;
                        h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * sin(phi(k, 1)) * cos(theta(k, 1)) + (ny-0.5*(nel(n_ar_in, 1)-1)) * sin(theta(k, 1))));
                    end
                end
                switch n_ar_in
                    case 1
                        H1(:, kx) = h * beta(k, 1);
                    case 2
                        H2(:, kx) = h * beta(k, 1);
                    case 3
                        H3(:, kx) = h * beta(k, 1);
                end
            end
        end
        user_number_group = floor(K / pilot_n);%5
        rp1 = zeros(n_arr(1, 1), pilot_n);
        rp2 = zeros(n_arr(2, 1), pilot_n);
        rp3 = zeros(n_arr(3, 1), pilot_n);
        for user_group_in = 1 : pilot_n%1:2
            for nn = 1 : user_number_group%5
                in_temp = pilot_n * (nn-1) + user_group_in;
                rp1(:, user_group_in) = rp1(:, user_group_in) + H1(:, in_temp);
                rp2(:, user_group_in) = rp2(:, user_group_in) + H2(:, in_temp);
                rp3(:, user_group_in) = rp3(:, user_group_in) + H3(:, in_temp);
            end
            if user_group_in == pilot_n && user_number_group * pilot_n < K
                for kk_temp = user_number_group * pilot_n + 1 : K
                    rp1(:, user_group_in) = rp1(:, user_group_in) + H1(:, kk_temp);
                    rp2(:, user_group_in) = rp2(:, user_group_in) + H2(:, kk_temp);
                    rp3(:, user_group_in) = rp3(:, user_group_in) + H3(:, kk_temp);
                end
            else
            end
            rp1(:, user_group_in) = (U1' * (sqrt(rho) * rp1(:, user_group_in) + randn(n_arr(1, 1), 1)));
            rp2(:, user_group_in) = (U2' * (sqrt(rho) * rp2(:, user_group_in) + randn(n_arr(2, 1), 1)));
            rp3(:, user_group_in) = (U3' * (sqrt(rho) * rp3(:, user_group_in) + randn(n_arr(3, 1), 1)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        corbes;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %generate beams, DOA estimation
        gbdoae;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %channel estimation, capacity calculation
        cecac;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %traditional processing
        tradipro;
        disp([signalpo_n, ii])
    end
    capacity(signalpo_n, 1:3) = capacity(signalpo_n, 1:3) * N_1;
end
capacity = capacity / Nite;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
capacity(:, 4:6) = capacity(:, 4:6) * N;
capacity = capacity / Tc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
set(h,'PaperType','A4');
xx = axes('FontSize',16);
plot(pilot_s, capacity(1:length(pilot_s), 1), 'b-d','LineWidth',2,'MarkerSize',14)
hold on
plot(pilot_s, capacity(1:length(pilot_s), 4), 'b--x','LineWidth',2,'MarkerSize',10)
plot(pilot_s, capacity(1:length(pilot_s), 2), 'r-d','LineWidth',2,'MarkerSize',14)
plot(pilot_s, capacity(1:length(pilot_s), 5), 'r--x','LineWidth',2,'MarkerSize',10)
plot(pilot_s, capacity(1:length(pilot_s), 3), 'k-d','LineWidth',2,'MarkerSize',14)
plot(pilot_s, capacity(1:length(pilot_s), 6), 'k--x','LineWidth',2,'MarkerSize',10)
ylim([min(capacity(:,1))-5, max(capacity(:,3))+30])
xlim([min(pilot_s), max(pilot_s)])
le = legend('proposed, N=252','traditional, N=252','proposed, N=1216','traditional, N=1216','proposed, N=4832','traditional, N=4832', 'Location', 'northwest');
set(le,'Fontsize',16,'Fontname','Times')
set(gca,'XTick',pilot_s)
xlabel('Number of the pilot symbols','Fontsize',20,'Fontname','Times')
ylabel('Sum rate (bps/Hz)','Fontsize',20,'Fontname','Times')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(h,'-dpdf','capacity_pilot')