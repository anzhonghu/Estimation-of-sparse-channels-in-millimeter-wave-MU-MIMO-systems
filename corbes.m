%correaltion



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%beam selection
temp_index1 = zeros(n_arr(1, 1), 4);
temp_index1(:, 3) = ones(n_arr(1, 1), 1);
sum_temp = 0;
for user_group_in = 1 : pilot_n%1:2
    if user_group_in < pilot_n
        eff_num = user_number_group;
    else
        eff_num = K - user_number_group * (pilot_n-1);
    end
    for nx = 1 : naz(1, 1)
        for ny = 1 : nel(1, 1)
            n = (ny - 1) * naz(1, 1) + nx;
            Cor1(nx, ny) = abs(rp1(n, user_group_in));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%111
    for nx = 2 : naz(1, 1)-1
        for ny = 2 : nel(1, 1)-1
            n = (ny - 1) * naz(1, 1) + nx;
            if n>=in1_l && n<=in1_u && 0 ~= temp_index1(n, 3)
                if Cor1(nx, ny)>Cor1(nx-1, ny) && Cor1(nx, ny)>Cor1(nx-1, ny-1) && Cor1(nx, ny)>Cor1(nx-1, ny+1) && Cor1(nx, ny)>Cor1(nx, ny-1)...
                        && Cor1(nx, ny)>Cor1(nx, ny+1) && Cor1(nx, ny)>Cor1(nx+1, ny) && Cor1(nx, ny)>Cor1(nx+1, ny-1) && Cor1(nx, ny)>Cor1(nx+1, ny+1)
                    temp_index1(n, 1) = 1;
                    temp_index1(n, 2) = Cor1(nx, ny);
                    temp_index1(n, 3) = 0;
                    temp_index1(n, 4) = temp_index1(n, 2);
                else
                    temp_index1(n, 3) = Cor1(nx, ny);
                end
            else
                temp_index1(n, 3) = 0;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%222
    nx = 1;
    for ny = 2 : nel(1, 1)-1
        n = (ny - 1) * naz(1, 1) + nx;
        if n>=in1_l && n<=in1_u && 0 ~= temp_index1(n, 3)
            if  Cor1(nx, ny)>Cor1(nx, ny-1)...
                    && Cor1(nx, ny)>Cor1(nx, ny+1) && Cor1(nx, ny)>Cor1(nx+1, ny) && Cor1(nx, ny)>Cor1(nx+1, ny-1) && Cor1(nx, ny)>Cor1(nx+1, ny+1)
                temp_index1(n, 1) = 1;
                temp_index1(n, 2) = Cor1(nx, ny);
                temp_index1(n, 3) = 0;
                temp_index1(n, 4) = temp_index1(n, 2);
            else
                temp_index1(n, 3) = Cor1(nx, ny);
            end
        else
            temp_index1(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%3333
    nx = naz(1, 1);
    for ny = 2 : nel(1, 1)-1
        n = (ny - 1) * naz(1, 1) + nx;
        if n>=in1_l && n<=in1_u && 0 ~= temp_index1(n, 3)
            if Cor1(nx, ny)>Cor1(nx-1, ny) && Cor1(nx, ny)>Cor1(nx-1, ny-1) && Cor1(nx, ny)>Cor1(nx-1, ny+1) && Cor1(nx, ny)>Cor1(nx, ny-1)...
                    && Cor1(nx, ny)>Cor1(nx, ny+1)
                temp_index1(n, 1) = 1;
                temp_index1(n, 2) = Cor1(nx, ny);
                temp_index1(n, 3) = 0;
                temp_index1(n, 4) = temp_index1(n, 2);
            else
                temp_index1(n, 3) = Cor1(nx, ny);
            end
        else
            temp_index1(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4444
    ny = 1;
    for nx = 2 : naz(1, 1)-1
        n = (ny - 1) * naz(1, 1) + nx;
        if n>=in1_l && n<=in1_u && 0 ~= temp_index1(n, 3)
            if Cor1(nx, ny)>Cor1(nx-1, ny)  && Cor1(nx, ny)>Cor1(nx-1, ny+1) ...
                    && Cor1(nx, ny)>Cor1(nx, ny+1) && Cor1(nx, ny)>Cor1(nx+1, ny)  && Cor1(nx, ny)>Cor1(nx+1, ny+1)
                temp_index1(n, 1) = 1;
                temp_index1(n, 2) = Cor1(nx, ny);
                temp_index1(n, 3) = 0;
                temp_index1(n, 4) = temp_index1(n, 2);
            else
                temp_index1(n, 3) = Cor1(nx, ny);
            end
        else
            temp_index1(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555555555555555
    ny = nel(1, 1);
    for nx = 2 : naz(1, 1)-1
        n = (ny - 1) * naz(1, 1) + nx;
        if n>=in1_l && n<=in1_u && 0 ~= temp_index1(n, 3)
            if Cor1(nx, ny)>Cor1(nx-1, ny) && Cor1(nx, ny)>Cor1(nx-1, ny-1) && Cor1(nx, ny)>Cor1(nx, ny-1)...
                    &&  Cor1(nx, ny)>Cor1(nx+1, ny) && Cor1(nx, ny)>Cor1(nx+1, ny-1)
                temp_index1(n, 1) = 1;
                temp_index1(n, 2) = Cor1(nx, ny);
                temp_index1(n, 3) = 0;
                temp_index1(n, 4) = temp_index1(n, 2);
            else
                temp_index1(n, 3) = Cor1(nx, ny);
            end
        else
            temp_index1(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(temp_index1(:, 1)) - sum_temp < eff_num
        [temp_value, temp_index_s] = sort(temp_index1(:, 3), 'descend');
        for n_com = 1 : eff_num - sum(temp_index1(:, 1)) + sum_temp
            temp_index1(temp_index_s(n_com, 1), 1) = 1;
            temp_index1(temp_index_s(n_com, 1), 2) = temp_value(n_com, 1);
            temp_index1(temp_index_s(n_com, 1), 3) = 0;
        end
    else
        if sum(temp_index1(:, 1)) - sum_temp > eff_num
            [temp_value_a, temp_index_sa] = sort(temp_index1(:, 4), 'descend');
            for n_com = eff_num+1 : sum(temp_index1(:, 1)) - sum_temp
                temp_index1(temp_index_sa(n_com, 1), 1) = 0;
                temp_index1(temp_index_sa(n_com, 1), 3) = 1;
            end
            temp_index1(:, 4) = 0;
        else
        end
    end
    sum_temp = sum(temp_index1(:, 1));
end
C1 = zeros(Mb, 2);
nk = 0;
index_label1 = zeros(n_arr(1, 1), 2);
for nx = 1 : naz(1, 1)
    for ny = 1 : nel(1, 1)
        n = (ny - 1) * naz(1, 1) + nx;
        if 0 == temp_index1(n, 1)
        else
            nk = nk + 1;
            C1(nk, 1) = n;
            C1(nk, 2) = temp_index1(n, 2);
            index_label1(n, 1) = 1;
        end
    end
end
[~, C_index] = sort(C1(1:K,1));
C1(1:K, :) = C1(C_index, :);
for k = 1 : K
    if abs(C1(k, 1)/naz(1, 1) - floor(C1(k, 1)/naz(1, 1)))<1e-2
        x_temp = naz(1, 1);
        y_temp = floor(C1(k, 1)/naz(1, 1));
    else
        y_temp = floor(C1(k, 1)/naz(1, 1)) + 1;
        x_temp = C1(k, 1) - (y_temp-1) * naz(1, 1);
    end
    for nx = 1 : naz(1, 1)
        for ny = 1 : nel(1, 1)
            n = (ny - 1) * naz(1, 1) + nx;
            if 0 == index_label1(n, 1)
                index_label1(n, 2) = abs(nx - x_temp) + abs(ny - y_temp);
            else
                index_label1(n, 2) = n_arr(1, 1);
            end
        end
    end
    [~, index_sortan] = sort(index_label1(:, 2));
    C1(K+(k-1)*(beam_numberforMS-1)+1:K+k*(beam_numberforMS-1), 1) = index_sortan(1:(beam_numberforMS-1), 1);
    C1(K+(k-1)*(beam_numberforMS-1)+1:K+k*(beam_numberforMS-1), 2) = temp_index1(index_sortan(1:(beam_numberforMS-1), 1), 2);
    index_label1(index_sortan(1:(beam_numberforMS-1), 1), 1) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_index2 = zeros(n_arr(2, 1), 4);
temp_index2(:, 3) = ones(n_arr(2, 1), 1);
sum_temp = 0;
for user_group_in = 1 : pilot_n%1:2
    if user_group_in < pilot_n
        eff_num = user_number_group;
    else
        eff_num = K - user_number_group * (pilot_n-1);
    end
    for nx = 1 : naz(2, 1)
        for ny = 1 : nel(2, 1)
            n = (ny - 1) * naz(2, 1) + nx;
            Cor2(nx, ny) = abs(rp2(n, user_group_in));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%111111
    for nx = 2 : naz(2, 1)-1
        for ny = 2 : nel(2, 1)-1
            n = (ny - 1) * naz(2, 1) + nx;
            if n>=in2_l && n<=in2_u && 0 ~= temp_index2(n, 3)
                if Cor2(nx, ny)>Cor2(nx-1, ny) && Cor2(nx, ny)>Cor2(nx-1, ny-1) && Cor2(nx, ny)>Cor2(nx-1, ny+1) && Cor2(nx, ny)>Cor2(nx, ny-1)...
                        && Cor2(nx, ny)>Cor2(nx, ny+1) && Cor2(nx, ny)>Cor2(nx+1, ny) && Cor2(nx, ny)>Cor2(nx+1, ny-1) && Cor2(nx, ny)>Cor2(nx+1, ny+1)
                    temp_index2(n, 1) = 1;
                    temp_index2(n, 2) = Cor2(nx, ny);
                    temp_index2(n, 3) = 0;
                    temp_index2(n, 4) = temp_index2(n, 2);
                else
                    temp_index2(n, 3) = Cor2(nx, ny);
                end
            else
                temp_index2(n, 3) = 0;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%222222222222222
    nx = 1;
    for ny = 2 : nel(2, 1)-1
        n = (ny - 1) * naz(2, 1) + nx;
        if n>=in2_l && n<=in2_u && 0 ~= temp_index2(n, 3)
            if  Cor2(nx, ny)>Cor2(nx, ny-1)...
                    && Cor2(nx, ny)>Cor2(nx, ny+1) && Cor2(nx, ny)>Cor2(nx+1, ny) && Cor2(nx, ny)>Cor2(nx+1, ny-1) && Cor2(nx, ny)>Cor2(nx+1, ny+1)
                temp_index2(n, 1) = 1;
                temp_index2(n, 2) = Cor2(nx, ny);
                temp_index2(n, 3) = 0;
                temp_index2(n, 4) = temp_index2(n, 2);
            else
                temp_index2(n, 3) = Cor2(nx, ny);
            end
        else
            temp_index2(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%333333
    nx = naz(2, 1);
    for ny = 2 : nel(2, 1)-1
        n = (ny - 1) * naz(2, 1) + nx;
        if n>=in2_l && n<=in2_u && 0 ~= temp_index2(n, 3)
            if Cor2(nx, ny)>Cor2(nx-1, ny) && Cor2(nx, ny)>Cor2(nx-1, ny-1) && Cor2(nx, ny)>Cor2(nx-1, ny+1) && Cor2(nx, ny)>Cor2(nx, ny-1)...
                    && Cor2(nx, ny)>Cor2(nx, ny+1)
                temp_index2(n, 1) = 1;
                temp_index2(n, 2) = Cor2(nx, ny);
                temp_index2(n, 3) = 0;
                temp_index2(n, 4) = temp_index2(n, 2);
            else
                temp_index2(n, 3) = Cor2(nx, ny);
            end
        else
            temp_index2(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%4444444444
    ny = 1;
    for nx = 2 : naz(2, 1)-1
        n = (ny - 1) * naz(2, 1) + nx;
        if n>=in2_l && n<=in2_u && 0 ~= temp_index2(n, 3)
            if Cor2(nx, ny)>Cor2(nx-1, ny)  && Cor2(nx, ny)>Cor2(nx-1, ny+1) ...
                    && Cor2(nx, ny)>Cor2(nx, ny+1) && Cor2(nx, ny)>Cor2(nx+1, ny)  && Cor2(nx, ny)>Cor2(nx+1, ny+1)
                temp_index2(n, 1) = 1;
                temp_index2(n, 2) = Cor2(nx, ny);
                temp_index2(n, 3) = 0;
                temp_index2(n, 4) = temp_index2(n, 2);
            else
                temp_index2(n, 3) = Cor2(nx, ny);
            end
        else
            temp_index2(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%555555
    ny = nel(2, 1);
    for nx = 2 : naz(2, 1)-1
        n = (ny - 1) * naz(2, 1) + nx;
        if n>=in2_l && n<=in2_u && 0 ~= temp_index2(n, 3)
            if Cor2(nx, ny)>Cor2(nx-1, ny) && Cor2(nx, ny)>Cor2(nx-1, ny-1)  && Cor2(nx, ny)>Cor2(nx, ny-1)...
                    && Cor2(nx, ny)>Cor2(nx+1, ny) && Cor2(nx, ny)>Cor2(nx+1, ny-1)
                temp_index2(n, 1) = 1;
                temp_index2(n, 2) = Cor2(nx, ny);
                temp_index2(n, 3) = 0;
                temp_index2(n, 4) = temp_index2(n, 2);
            else
                temp_index2(n, 3) = Cor2(nx, ny);
            end
        else
            temp_index2(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(temp_index2(:, 1)) - sum_temp < eff_num
        [temp_value, temp_index_s] = sort(temp_index2(:, 3), 'descend');
        for n_com = 1 : eff_num-sum(temp_index2(:, 1))+sum_temp
            temp_index2(temp_index_s(n_com, 1), 1) = 1;
            temp_index2(temp_index_s(n_com, 1), 2) = temp_value(n_com, 1);
            temp_index2(temp_index_s(n_com, 1), 3) = 0;
        end
    else
        if sum(temp_index2(:, 1)) - sum_temp > eff_num
            [temp_value_a, temp_index_sa] = sort(temp_index2(:, 4), 'descend');
            for n_com = eff_num+1 : sum(temp_index2(:, 1))-sum_temp
                temp_index2(temp_index_sa(n_com, 1), 1) = 0;
                temp_index2(temp_index_sa(n_com, 1), 3) = 1;
            end
            temp_index2(:, 4) = 0;
        else
        end
    end
    sum_temp = sum(temp_index2(:, 1));
end
C2 = zeros(Mb, 2);
nk = 0;
index_label2 = zeros(n_arr(2, 1), 2);
for nx = 1 : naz(2, 1)
    for ny = 1 : nel(2, 1)
        n = (ny - 1) * naz(2, 1) + nx;
        if 0 == temp_index2(n, 1)
        else
            nk = nk + 1;
            C2(nk, 1) = n;
            C2(nk, 2) = temp_index2(n, 2);
            index_label2(n, 1) = 1;
        end
    end
end
[~, C_index] = sort(C2(1:K,1));
C2(1:K, :) = C2(C_index, :);
for k = 1 : K
    if abs(C2(k, 1)/naz(2, 1) - floor(C2(k, 1)/naz(2, 1)))<1e-2
        x_temp = naz(2, 1);
        y_temp = floor(C2(k, 1)/naz(2, 1));
    else
        y_temp = floor(C2(k, 1)/naz(2, 1)) + 1;
        x_temp = C2(k, 1) - (y_temp-1) * naz(2, 1);
    end
    for nx = 1 : naz(2, 1)
        for ny = 1 : nel(2, 1)
            n = (ny - 1) * naz(2, 1) + nx;
            if 0 == index_label2(n, 1)
                index_label2(n, 2) = abs(nx - x_temp) + abs(ny - y_temp);
            else
                index_label2(n, 2) = n_arr(2, 1);
            end
        end
    end
    [~, index_sortan] = sort(index_label2(:, 2));
    C2(K+(k-1)*(beam_numberforMS-1)+1:K+k*(beam_numberforMS-1), 1) = index_sortan(1:(beam_numberforMS-1), 1);
    C2(K+(k-1)*(beam_numberforMS-1)+1:K+k*(beam_numberforMS-1), 2) = temp_index2(index_sortan(1:(beam_numberforMS-1), 1), 2);
    index_label2(index_sortan(1:(beam_numberforMS-1), 1), 1) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_index3 = zeros(n_arr(3, 1), 4);
temp_index3(:, 3) = ones(n_arr(3, 1), 1);
sum_temp = 0;
for user_group_in = 1 : pilot_n%1:2
    if user_group_in < pilot_n
        eff_num = user_number_group;
    else
        eff_num = K - user_number_group * (pilot_n-1);
    end
    for nx = 1 : naz(3, 1)
        for ny = 1 : nel(3, 1)
            n = (ny - 1) * naz(3, 1) + nx;
            Cor3(nx, ny) = abs(rp3(n, user_group_in));
        end
    end
    %%%%%%%%%%%%%%%%%%%1111111
    for nx = 2 : naz(3, 1)-1
        for ny = 2 : nel(3, 1)-1
            n = (ny - 1) * naz(3, 1) + nx;
            if n>=in3_l && n<=in3_u && 0 ~= temp_index3(n, 3)
                if Cor3(nx, ny)>Cor3(nx-1, ny) && Cor3(nx, ny)>Cor3(nx-1, ny-1) && Cor3(nx, ny)>Cor3(nx-1, ny+1) && Cor3(nx, ny)>Cor3(nx, ny-1)...
                        && Cor3(nx, ny)>Cor3(nx, ny+1) && Cor3(nx, ny)>Cor3(nx+1, ny) && Cor3(nx, ny)>Cor3(nx+1, ny-1) && Cor3(nx, ny)>Cor3(nx+1, ny+1)
                    temp_index3(n, 1) = 1;
                    temp_index3(n, 2) = Cor3(nx, ny);
                    temp_index3(n, 3) = 0;
                    temp_index3(n, 4) = temp_index3(n, 2);
                else
                    temp_index3(n, 3) = Cor3(nx, ny);
                end
            else
                temp_index3(n, 3) = 0;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%22222222
    nx = 1;
    for ny = 2 : nel(3, 1)-1
        n = (ny - 1) * naz(3, 1) + nx;
        if n>=in3_l && n<=in3_u && 0 ~= temp_index3(n, 3)
            if  Cor3(nx, ny)>Cor3(nx, ny-1)...
                    && Cor3(nx, ny)>Cor3(nx, ny+1) && Cor3(nx, ny)>Cor3(nx+1, ny) && Cor3(nx, ny)>Cor3(nx+1, ny-1) && Cor3(nx, ny)>Cor3(nx+1, ny+1)
                temp_index3(n, 1) = 1;
                temp_index3(n, 2) = Cor3(nx, ny);
                temp_index3(n, 3) = 0;
                temp_index3(n, 4) = temp_index3(n, 2);
            else
                temp_index3(n, 3) = Cor3(nx, ny);
            end
        else
            temp_index3(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%33333
    nx = naz(3, 1);
    for ny = 2 : nel(3, 1)-1
        n = (ny - 1) * naz(3, 1) + nx;
        if n>=in3_l && n<=in3_u && 0 ~= temp_index3(n, 3)
            if Cor3(nx, ny)>Cor3(nx-1, ny) && Cor3(nx, ny)>Cor3(nx-1, ny-1) && Cor3(nx, ny)>Cor3(nx-1, ny+1) && Cor3(nx, ny)>Cor3(nx, ny-1)...
                    && Cor3(nx, ny)>Cor3(nx, ny+1)
                temp_index3(n, 1) = 1;
                temp_index3(n, 2) = Cor3(nx, ny);
                temp_index3(n, 3) = 0;
                temp_index3(n, 4) = temp_index3(n, 2);
            else
                temp_index3(n, 3) = Cor3(nx, ny);
            end
        else
            temp_index3(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%4444444
    ny = 1;
    for nx = 2 : naz(3, 1)-1
        n = (ny - 1) * naz(3, 1) + nx;
        if n>=in3_l && n<=in3_u && 0 ~= temp_index3(n, 3)
            if Cor3(nx, ny)>Cor3(nx-1, ny) && Cor3(nx, ny)>Cor3(nx-1, ny+1) ...
                    && Cor3(nx, ny)>Cor3(nx, ny+1) && Cor3(nx, ny)>Cor3(nx+1, ny) && Cor3(nx, ny)>Cor3(nx+1, ny+1)
                temp_index3(n, 1) = 1;
                temp_index3(n, 2) = Cor3(nx, ny);
                temp_index3(n, 3) = 0;
                temp_index3(n, 4) = temp_index3(n, 2);
            else
                temp_index3(n, 3) = Cor3(nx, ny);
            end
        else
            temp_index3(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%55555555555
    ny = nel(3, 1);
    for nx = 2 : naz(3, 1)-1
        n = (ny - 1) * naz(3, 1) + nx;
        if n>=in3_l && n<=in3_u && 0 ~= temp_index3(n, 3)
            if Cor3(nx, ny)>Cor3(nx-1, ny) && Cor3(nx, ny)>Cor3(nx-1, ny-1) && Cor3(nx, ny)>Cor3(nx, ny-1)...
                    && Cor3(nx, ny)>Cor3(nx+1, ny) && Cor3(nx, ny)>Cor3(nx+1, ny-1)
                temp_index3(n, 1) = 1;
                temp_index3(n, 2) = Cor3(nx, ny);
                temp_index3(n, 3) = 0;
                temp_index3(n, 4) = temp_index3(n, 2);
            else
                temp_index3(n, 3) = Cor3(nx, ny);
            end
        else
            temp_index3(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(temp_index3(:, 1)) - sum_temp < eff_num
        [temp_value, temp_index_s] = sort(temp_index3(:, 3), 'descend');
        for n_com = 1 : eff_num-sum(temp_index3(:, 1))+sum_temp
            temp_index3(temp_index_s(n_com, 1), 1) = 1;
            temp_index3(temp_index_s(n_com, 1), 2) = temp_value(n_com, 1);
            temp_index3(temp_index_s(n_com, 1), 3) = 0;
        end
    else
        if sum(temp_index3(:, 1)) - sum_temp > eff_num
            [temp_value_a, temp_index_sa] = sort(temp_index3(:, 4), 'descend');
            for n_com = eff_num+1 : sum(temp_index3(:, 1))-sum_temp
                temp_index3(temp_index_sa(n_com, 1), 1) = 0;
                temp_index3(temp_index_sa(n_com, 1), 3) = 1;
            end
          temp_index3(:, 4) = 0;
        else
        end
    end
    sum_temp = sum(temp_index3(:, 1));
end
C3 = zeros(Mb, 2);
nk = 0;
index_label3 = zeros(n_arr(3, 1), 2);
for nx = 1 : naz(3, 1)
    for ny = 1 : nel(3, 1)
        n = (ny - 1) * naz(3, 1) + nx;
        if 0 == temp_index3(n, 1)
        else
            nk = nk + 1;
            C3(nk, 1) = n;
            C3(nk, 2) = temp_index3(n, 2);
            index_label3(n, 1) = 1;
        end
    end
end
[~, C_index] = sort(C3(1:K,1));
C3(1:K, :) = C3(C_index, :);
for k = 1 : K
    if abs(C3(k, 1)/naz(3, 1) - floor(C3(k, 1)/naz(3, 1)))<1e-2
        x_temp = naz(3, 1);
        y_temp = floor(C3(k, 1)/naz(3, 1));
    else
        y_temp = floor(C3(k, 1)/naz(3, 1)) + 1;
        x_temp = C3(k, 1) - (y_temp-1) * naz(3, 1);
    end
    for nx = 1 : naz(3, 1)
        for ny = 1 : nel(3, 1)
            n = (ny - 1) * naz(3, 1) + nx;
            if 0 == index_label3(n, 1)
                index_label3(n, 2) = abs(nx - x_temp) + abs(ny - y_temp);
            else
                index_label3(n, 2) = n_arr(3, 1);
            end
        end
    end
    [~, index_sortan] = sort(index_label3(:, 2));
    C3(K+(k-1)*(beam_numberforMS-1)+1:K+k*(beam_numberforMS-1), 1) = index_sortan(1:(beam_numberforMS-1), 1);
    C3(K+(k-1)*(beam_numberforMS-1)+1:K+k*(beam_numberforMS-1), 2) = temp_index3(index_sortan(1:(beam_numberforMS-1), 1), 2);
    index_label3(index_sortan(1:(beam_numberforMS-1), 1), 1) = 1;
end