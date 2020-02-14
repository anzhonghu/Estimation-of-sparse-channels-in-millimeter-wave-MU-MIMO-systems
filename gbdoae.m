        U1b = zeros(n_arr(1, 1), Mb);
        spatailf1 = zeros(K, 3);
        for k = 1 : K
            if k > user_number_group * pilot_n
                gr_n = pilot_n;
            else
                gr_n = ceil(k / user_number_group);
            end
            U1b(:, (k-1)*beam_numberforMS+1) = U1(:, C1(k, 1));
            flag = 0;
            spa_st = zeros(beam_numberforMS, 2);
            if abs(C1(k, 1)/naz(1, 1) - floor(C1(k, 1)/naz(1, 1)))<1e-2
                x_t = naz(1, 1);
                y_t = floor(C1(k, 1)/naz(1, 1));
            else
                y_t = floor(C1(k, 1)/naz(1, 1)) + 1;
                x_t = C1(k, 1) - (y_t-1) * naz(1, 1);
            end
            spa_st(1,1) = x_t;
            spa_st(1,2) = y_t;
            for nnnan = 1 : beam_numberforMS-1
                if abs(C1(K+(k-1)*(beam_numberforMS-1)+nnnan, 1)/naz(1, 1) - floor(C1(K+(k-1)*(beam_numberforMS-1)+nnnan, 1)/naz(1, 1)))<1e-2
                    x_temp = naz(1, 1);
                    y_temp = floor(C1(K+(k-1)*(beam_numberforMS-1)+nnnan, 1)/naz(1, 1));
                else
                    y_temp = floor(C1(K+(k-1)*(beam_numberforMS-1)+nnnan, 1)/naz(1, 1)) + 1;
                    x_temp = C1(K+(k-1)*(beam_numberforMS-1)+nnnan, 1) - (y_temp-1) * naz(1, 1);
                end
                U1b(:, (k-1)*beam_numberforMS+1+nnnan) = U1(:, C1(K+(k-1)*(beam_numberforMS-1)+nnnan, 1));
                spa_st(nnnan+1, :) = [x_temp, y_temp];
            end
            x_t = x_t/naz(1, 1)*2-1;
            y_t = y_t/nel(1, 1)*2-1;
            spa_st(:, 1) = spa_st(:, 1)/naz(1, 1)*2-1;
            spa_st(:, 2) = spa_st(:, 2)/nel(1, 1)*2-1;
            spa_diff = zeros(Spad*Spad, 1);
            for nnnan = 2 : beam_numberforMS
                if abs(x_t-spa_st(nnnan, 1))>2/naz(1, 1) || abs(y_t-spa_st(nnnan, 2))>2/nel(1, 1) || 1==flag
                else
                    flag = 1;
                    for spa_cx = 1 : Spad
                        for spa_cy = 1 : Spad
                            spa_x_temp = min(spa_st(:,1)) + (max(spa_st(:,1))-min(spa_st(:,1))) / Spad * spa_cx;
                            spa_y_temp = min(spa_st(:,2)) + (max(spa_st(:,2))-min(spa_st(:,2))) / Spad * spa_cy;
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if spa_st(nnnan, 1) == spa_x_temp && spa_st(nnnan, 2) == spa_y_temp
                                lambdamk = n_arr(1, 1);
                            else
                            end
                            if spa_st(nnnan, 1) ~= spa_x_temp && spa_st(nnnan, 2) == spa_y_temp
                                lambdamk = nel(1, 1) * (1-exp(1i*pi*naz(1, 1)*(spa_st(nnnan, 1)-spa_x_temp))) / (1-exp(1i*pi*(spa_st(nnnan, 1)-spa_x_temp)));
                            else
                            end
                            if spa_st(nnnan, 1) == spa_x_temp && spa_st(nnnan, 2) ~= spa_y_temp
                                lambdamk = naz(1, 1) * (1-exp(1i*pi*nel(1, 1)*(spa_st(nnnan, 2)-spa_y_temp))) / (1-exp(1i*pi*(spa_st(nnnan, 2)-spa_y_temp)));
                            else
                            end
                            if spa_st(nnnan, 1) ~= spa_x_temp && spa_st(nnnan, 2) ~= spa_y_temp
                                lambdamk = (1-exp(1i*pi*naz(1, 1)*(spa_st(nnnan, 1)-spa_x_temp))) * (1-exp(1i*pi*nel(1, 1)*(spa_st(nnnan, 2)-spa_y_temp)))...
                                / (1-exp(1i*pi*(spa_st(nnnan, 1)-spa_x_temp))) / (1-exp(1i*pi*(spa_st(nnnan, 2)-spa_y_temp)));
                            else
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if x_t == spa_x_temp && y_t == spa_y_temp
                                lambdackk = n_arr(1, 1);
                            else
                            end
                            if x_t ~= spa_x_temp && y_t == spa_y_temp
                                lambdackk = nel(1, 1) * (1-exp(1i*pi*naz(1, 1)*(x_t-spa_x_temp))) / (1-exp(1i*pi*(x_t-spa_x_temp)));
                            else
                            end
                            if x_t == spa_x_temp && y_t ~= spa_y_temp
                                lambdackk = naz(1, 1) * (1-exp(1i*pi*nel(1, 1)*(y_t-spa_y_temp))) / (1-exp(1i*pi*(y_t-spa_y_temp)));
                            else
                            end
                            if x_t ~= spa_x_temp && y_t ~= spa_y_temp
                                lambdackk = (1-exp(1i*pi*naz(1, 1)*(x_t-spa_x_temp))) * (1-exp(1i*pi*nel(1, 1)*(y_t-spa_y_temp)))...
                                / (1-exp(1i*pi*(x_t-spa_x_temp))) / (1-exp(1i*pi*(y_t-spa_y_temp)));
                            else
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            spa_diff((spa_cx-1)*Spad+spa_cy, 1) = abs(rp1(C1(K+(k-1)*(beam_numberforMS-1)+nnnan-1, 1), gr_n)/rp1(C1(k, 1), gr_n)-lambdamk/lambdackk);
                        end
                    end
                    [~,spa_diff_ind] = sort(spa_diff);
                    if abs(spa_diff_ind(1,1)/Spad - floor(spa_diff_ind(1,1)/Spad))<1e-2
                        spa_cy = Spad;
                        spa_cx = floor(spa_diff_ind(1,1)/Spad);
                    else
                        spa_cx = floor(spa_diff_ind(1,1)/Spad) + 1;
                        spa_cy = spa_diff_ind(1,1) - (spa_cx-1) * Spad;
                    end
                    spatailf1(k, 1) = min(spa_st(:,1)) + (max(spa_st(:,1))-min(spa_st(:,1))) / Spad * spa_cx;
                    spatailf1(k, 2) = min(spa_st(:,2)) + (max(spa_st(:,2))-min(spa_st(:,2))) / Spad * spa_cy;
                end
            end
            if 0 == flag
                spatailf1(k, 1) = x_t;
                spatailf1(k, 2) = y_t;
            else
            end
            lambdae = exp(1i*0.5*pi*((-naz(1, 1)+3)*spatailf1(k, 1)+(-nel(1, 1)+3)*spatailf1(k, 2)));
            if x_t == spatailf1(k, 1) && y_t == spatailf1(k, 2)
                lambdamk = n_arr(1, 1);
            else
            end
            if x_t ~= spatailf1(k, 1) && y_t == spatailf1(k, 2)
                lambdamk = nel(1, 1) * (1-exp(1i*pi*naz(1, 1)*(x_t-spatailf1(k, 1)))) / (1-exp(1i*pi*(x_t-spatailf1(k, 1))));
            else
            end
            if x_t == spatailf1(k, 1) && y_t ~= spatailf1(k, 2)
                lambdamk = naz(1, 1) * (1-exp(1i*pi*nel(1, 1)*(y_t-spatailf1(k, 2)))) / (1-exp(1i*pi*(y_t-spatailf1(k, 2))));
            else
            end
            if x_t ~= spatailf1(k, 1) && y_t ~= spatailf1(k, 2)
                lambdamk = (1-exp(1i*pi*naz(1, 1)*(x_t-spatailf1(k, 1)))) * (1-exp(1i*pi*nel(1, 1)*(y_t-spatailf1(k, 2))))...
                / (1-exp(1i*pi*(x_t-spatailf1(k, 1)))) / (1-exp(1i*pi*(y_t-spatailf1(k, 2))));
            else
            end
            spatailf1(k, 3) = sqrt(n_arr(1, 1)) * rp1(C1(k, 1), gr_n) / sqrt(rho) / lambdae / lambdamk;
        end
        U2b = zeros(n_arr(2, 1), Mb);
        spatailf2 = zeros(K, 3);
        for k = 1 : K
            U2b(:, (k-1)*beam_numberforMS+1) = U2(:, C2(k, 1));
            flag = 0;
            spa_st = zeros(beam_numberforMS, 2);
            if abs(C2(k, 1)/naz(2, 1) - floor(C2(k, 1)/naz(2, 1)))<1e-2
                x_t = naz(2, 1);
                y_t = floor(C2(k, 1)/naz(2, 1));
            else
                y_t = floor(C2(k, 1)/naz(2, 1)) + 1;
                x_t = C2(k, 1) - (y_t-1) * naz(2, 1);
            end
            spa_st(1,1) = x_t;
            spa_st(1,2) = y_t;
            for nnnan = 1 : beam_numberforMS-1
                if abs(C2(K+(k-1)*(beam_numberforMS-1)+nnnan, 1)/naz(2, 1) - floor(C2(K+(k-1)*(beam_numberforMS-1)+nnnan, 1)/naz(2, 1)))<1e-2
                    x_temp = naz(2, 1);
                    y_temp = floor(C2(K+(k-1)*(beam_numberforMS-1)+nnnan, 1)/naz(2, 1));
                else
                    y_temp = floor(C2(K+(k-1)*(beam_numberforMS-1)+nnnan, 1)/naz(2, 1)) + 1;
                    x_temp = C2(K+(k-1)*(beam_numberforMS-1)+nnnan, 1) - (y_temp-1) * naz(2, 1);
                end
                U2b(:, (k-1)*beam_numberforMS+1+nnnan) = U2(:, C2(K+(k-1)*(beam_numberforMS-1)+nnnan, 1));
                spa_st(nnnan+1, :) = [x_temp, y_temp];
            end
            x_t = x_t/naz(2, 1)*2-1;
            y_t = y_t/nel(2, 1)*2-1;
            spa_st(:, 1) = spa_st(:, 1)/naz(2, 1)*2-1;
            spa_st(:, 2) = spa_st(:, 2)/nel(2, 1)*2-1;
            spa_diff = zeros(Spad*Spad, 1);
            for nnnan = 2 : beam_numberforMS
                if abs(x_t-spa_st(nnnan, 1))>2/naz(2, 1) || abs(y_t-spa_st(nnnan, 2))>2/nel(2, 1) || 1==flag
                else
                    flag = 1;
                    for spa_cx = 1 : Spad
                        for spa_cy = 1 : Spad
                            spa_x_temp = min(spa_st(:,1)) + (max(spa_st(:,1))-min(spa_st(:,1))) / Spad * spa_cx;
                            spa_y_temp = min(spa_st(:,2)) + (max(spa_st(:,2))-min(spa_st(:,2))) / Spad * spa_cy;
                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if spa_st(nnnan, 1) == spa_x_temp && spa_st(nnnan, 2) == spa_y_temp
                                lambdamk = n_arr(2, 1);
                            else
                            end
                            if spa_st(nnnan, 1) ~= spa_x_temp && spa_st(nnnan, 2) == spa_y_temp
                                lambdamk = nel(2, 1) * (1-exp(1i*pi*naz(1, 1)*(spa_st(nnnan, 1)-spa_x_temp))) / (1-exp(1i*pi*(spa_st(nnnan, 1)-spa_x_temp)));
                            else
                            end
                            if spa_st(nnnan, 1) == spa_x_temp && spa_st(nnnan, 2) ~= spa_y_temp
                                lambdamk = naz(2, 1) * (1-exp(1i*pi*nel(1, 1)*(spa_st(nnnan, 2)-spa_y_temp))) / (1-exp(1i*pi*(spa_st(nnnan, 2)-spa_y_temp)));
                            else
                            end
                            if spa_st(nnnan, 1) ~= spa_x_temp && spa_st(nnnan, 2) ~= spa_y_temp
                                lambdamk = (1-exp(1i*pi*naz(1, 1)*(spa_st(nnnan, 1)-spa_x_temp))) * (1-exp(1i*pi*nel(1, 1)*(spa_st(nnnan, 2)-spa_y_temp)))...
                                / (1-exp(1i*pi*(spa_st(nnnan, 1)-spa_x_temp))) / (1-exp(1i*pi*(spa_st(nnnan, 2)-spa_y_temp)));
                            else
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if x_t == spa_x_temp && y_t == spa_y_temp
                                lambdackk = n_arr(2, 1);
                            else
                            end
                            if x_t ~= spa_x_temp && y_t == spa_y_temp
                                lambdackk = nel(2, 1) * (1-exp(1i*pi*naz(1, 1)*(x_t-spa_x_temp))) / (1-exp(1i*pi*(x_t-spa_x_temp)));
                            else
                            end
                            if x_t == spa_x_temp && y_t ~= spa_y_temp
                                lambdackk = naz(2, 1) * (1-exp(1i*pi*nel(1, 1)*(y_t-spa_y_temp))) / (1-exp(1i*pi*(y_t-spa_y_temp)));
                            else
                            end
                            if x_t ~= spa_x_temp && y_t ~= spa_y_temp
                                lambdackk = (1-exp(1i*pi*naz(1, 1)*(x_t-spa_x_temp))) * (1-exp(1i*pi*nel(1, 1)*(y_t-spa_y_temp)))...
                                / (1-exp(1i*pi*(x_t-spa_x_temp))) / (1-exp(1i*pi*(y_t-spa_y_temp)));
                            else
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            spa_diff((spa_cx-1)*Spad+spa_cy, 1) = abs(rp2(C2(K+(k-1)*(beam_numberforMS-1)+nnnan-1, 1), gr_n)/rp2(C2(k, 1), gr_n)-lambdamk/lambdackk);
                        end
                    end
                    [~,spa_diff_ind] = sort(spa_diff);
                    if abs(spa_diff_ind(1,1)/Spad - floor(spa_diff_ind(1,1)/Spad))<1e-2
                        spa_cy = Spad;
                        spa_cx = floor(spa_diff_ind(1,1)/Spad);
                    else
                        spa_cx = floor(spa_diff_ind(1,1)/Spad) + 1;
                        spa_cy = spa_diff_ind(1,1) - (spa_cx-1) * Spad;
                    end
                    spatailf2(k, 1) = min(spa_st(:,1)) + (max(spa_st(:,1))-min(spa_st(:,1))) / Spad * spa_cx;
                    spatailf2(k, 2) = min(spa_st(:,2)) + (max(spa_st(:,2))-min(spa_st(:,2))) / Spad * spa_cy;
                end
            end
            if 0 == flag
                spatailf2(k, 1) = x_t;
                spatailf2(k, 2) = y_t;
            else
            end
            lambdae = exp(1i*0.5*pi*((-naz(2, 1)+3)*spatailf2(k, 1)+(-nel(2, 1)+3)*spatailf2(k, 2)));
            if x_t == spatailf2(k, 1) && y_t == spatailf2(k, 2)
                lambdamk = n_arr(2, 1);
            else
            end
            if x_t ~= spatailf2(k, 1) && y_t == spatailf2(k, 2)
                lambdamk = nel(2, 1) * (1-exp(1i*pi*naz(2, 1)*(x_t-spatailf2(k, 1)))) / (1-exp(1i*pi*(x_t-spatailf2(k, 1))));
            else
            end
            if x_t == spatailf2(k, 1) && y_t ~= spatailf2(k, 2)
                lambdamk = naz(2, 1) * (1-exp(1i*pi*nel(2, 1)*(y_t-spatailf2(k, 2)))) / (1-exp(1i*pi*(y_t-spatailf2(k, 2))));
            else
            end
            if x_t ~= spatailf2(k, 1) && y_t ~= spatailf2(k, 2)
                lambdamk = (1-exp(1i*pi*naz(2, 1)*(x_t-spatailf2(k, 1)))) * (1-exp(1i*pi*nel(2, 1)*(y_t-spatailf2(k, 2))))...
                / (1-exp(1i*pi*(x_t-spatailf2(k, 1)))) / (1-exp(1i*pi*(y_t-spatailf2(k, 2))));
            else
            end
            spatailf2(k, 3) = sqrt(n_arr(2, 1)) * rp2(C2(k, 1), gr_n) / sqrt(rho) / lambdae / lambdamk;
        end
       U3b = zeros(n_arr(3, 1), Mb);
        spatailf3 = zeros(K, 3);
        for k = 1 : K
            U3b(:, (k-1)*beam_numberforMS+1) = U3(:, C3(k, 1));
            flag = 0;
            spa_st = zeros(beam_numberforMS, 2);
            if abs(C3(k, 1)/naz(3, 1) - floor(C3(k, 1)/naz(3, 1)))<1e-2
                x_t = naz(3, 1);
                y_t = floor(C3(k, 1)/naz(3, 1));
            else
                y_t = floor(C3(k, 1)/naz(3, 1)) + 1;
                x_t = C3(k, 1) - (y_t-1) * naz(3, 1);
            end
            spa_st(1,1) = x_t;
            spa_st(1,2) = y_t;
            for nnnan = 1 : beam_numberforMS-1
                if abs(C3(K+(k-1)*(beam_numberforMS-1)+nnnan, 1)/naz(3, 1) - floor(C3(K+(k-1)*(beam_numberforMS-1)+nnnan, 1)/naz(3, 1)))<1e-2
                    x_temp = naz(3, 1);
                    y_temp = floor(C3(K+(k-1)*(beam_numberforMS-1)+nnnan, 1)/naz(3, 1));
                else
                    y_temp = floor(C3(K+(k-1)*(beam_numberforMS-1)+nnnan, 1)/naz(3, 1)) + 1;
                    x_temp = C3(K+(k-1)*(beam_numberforMS-1)+nnnan, 1) - (y_temp-1) * naz(3, 1);
                end
                U3b(:, (k-1)*beam_numberforMS+1+nnnan) = U3(:, C3(K+(k-1)*(beam_numberforMS-1)+nnnan, 1));
                spa_st(nnnan+1, :) = [x_temp, y_temp];
            end
            x_t = x_t/naz(3, 1)*2-1;
            y_t = y_t/nel(3, 1)*2-1;
            spa_st(:, 1) = spa_st(:, 1)/naz(3, 1)*2-1;
            spa_st(:, 2) = spa_st(:, 2)/nel(3, 1)*2-1;
            spa_diff = zeros(Spad*Spad, 1);
            for nnnan = 2 : beam_numberforMS
                if abs(x_t-spa_st(nnnan, 1))>2/naz(3, 1) || abs(y_t-spa_st(nnnan, 2))>2/nel(3, 1) || 1==flag
                else
                    flag = 1;
                    for spa_cx = 1 : Spad
                        for spa_cy = 1 : Spad
                            spa_x_temp = min(spa_st(:,1)) + (max(spa_st(:,1))-min(spa_st(:,1))) / Spad * spa_cx;
                            spa_y_temp = min(spa_st(:,2)) + (max(spa_st(:,2))-min(spa_st(:,2))) / Spad * spa_cy;
                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if spa_st(nnnan, 1) == spa_x_temp && spa_st(nnnan, 2) == spa_y_temp
                                lambdamk = n_arr(3, 1);
                            else
                            end
                            if spa_st(nnnan, 1) ~= spa_x_temp && spa_st(nnnan, 2) == spa_y_temp
                                lambdamk = nel(3, 1) * (1-exp(1i*pi*naz(1, 1)*(spa_st(nnnan, 1)-spa_x_temp))) / (1-exp(1i*pi*(spa_st(nnnan, 1)-spa_x_temp)));
                            else
                            end
                            if spa_st(nnnan, 1) == spa_x_temp && spa_st(nnnan, 2) ~= spa_y_temp
                                lambdamk = naz(3, 1) * (1-exp(1i*pi*nel(1, 1)*(spa_st(nnnan, 2)-spa_y_temp))) / (1-exp(1i*pi*(spa_st(nnnan, 2)-spa_y_temp)));
                            else
                            end
                            if spa_st(nnnan, 1) ~= spa_x_temp && spa_st(nnnan, 2) ~= spa_y_temp
                                lambdamk = (1-exp(1i*pi*naz(1, 1)*(spa_st(nnnan, 1)-spa_x_temp))) * (1-exp(1i*pi*nel(1, 1)*(spa_st(nnnan, 2)-spa_y_temp)))...
                                / (1-exp(1i*pi*(spa_st(nnnan, 1)-spa_x_temp))) / (1-exp(1i*pi*(spa_st(nnnan, 2)-spa_y_temp)));
                            else
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if x_t == spa_x_temp && y_t == spa_y_temp
                                lambdackk = n_arr(3, 1);
                            else
                            end
                            if x_t ~= spa_x_temp && y_t == spa_y_temp
                                lambdackk = nel(3, 1) * (1-exp(1i*pi*naz(1, 1)*(x_t-spa_x_temp))) / (1-exp(1i*pi*(x_t-spa_x_temp)));
                            else
                            end
                            if x_t == spa_x_temp && y_t ~= spa_y_temp
                                lambdackk = naz(3, 1) * (1-exp(1i*pi*nel(1, 1)*(y_t-spa_y_temp))) / (1-exp(1i*pi*(y_t-spa_y_temp)));
                            else
                            end
                            if x_t ~= spa_x_temp && y_t ~= spa_y_temp
                                lambdackk = (1-exp(1i*pi*naz(1, 1)*(x_t-spa_x_temp))) * (1-exp(1i*pi*nel(1, 1)*(y_t-spa_y_temp)))...
                                / (1-exp(1i*pi*(x_t-spa_x_temp))) / (1-exp(1i*pi*(y_t-spa_y_temp)));
                            else
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            spa_diff((spa_cx-1)*Spad+spa_cy, 1) = abs(rp3(C3(K+(k-1)*(beam_numberforMS-1)+nnnan-1, 1), gr_n)/rp3(C3(k, 1), gr_n)-lambdamk/lambdackk);
                        end
                    end
                    [~,spa_diff_ind] = sort(spa_diff);
                    if abs(spa_diff_ind(1,1)/Spad - floor(spa_diff_ind(1,1)/Spad))<1e-2
                        spa_cy = Spad;
                        spa_cx = floor(spa_diff_ind(1,1)/Spad);
                    else
                        spa_cx = floor(spa_diff_ind(1,1)/Spad) + 1;
                        spa_cy = spa_diff_ind(1,1) - (spa_cx-1) * Spad;
                    end
                    spatailf3(k, 1) = min(spa_st(:,1)) + (max(spa_st(:,1))-min(spa_st(:,1))) / Spad * spa_cx;
                    spatailf3(k, 2) = min(spa_st(:,2)) + (max(spa_st(:,2))-min(spa_st(:,2))) / Spad * spa_cy;
                end
            end
            if 0 == flag
                spatailf3(k, 1) = x_t;
                spatailf3(k, 2) = y_t;
            else
            end
            lambdae = exp(1i*0.5*pi*((-naz(3, 1)+3)*spatailf3(k, 1)+(-nel(3, 1)+3)*spatailf3(k, 2)));
            if x_t == spatailf3(k, 1) && y_t == spatailf3(k, 2)
                lambdamk = n_arr(3, 1);
            else
            end
            if x_t ~= spatailf3(k, 1) && y_t == spatailf3(k, 2)
                lambdamk = nel(3, 1) * (1-exp(1i*pi*naz(3, 1)*(x_t-spatailf3(k, 1)))) / (1-exp(1i*pi*(x_t-spatailf3(k, 1))));
            else
            end
            if x_t == spatailf3(k, 1) && y_t ~= spatailf3(k, 2)
                lambdamk = naz(3, 1) * (1-exp(1i*pi*nel(3, 1)*(y_t-spatailf3(k, 2)))) / (1-exp(1i*pi*(y_t-spatailf3(k, 2))));
            else
            end
            if x_t ~= spatailf3(k, 1) && y_t ~= spatailf3(k, 2)
                lambdamk = (1-exp(1i*pi*naz(3, 1)*(x_t-spatailf3(k, 1)))) * (1-exp(1i*pi*nel(3, 1)*(y_t-spatailf3(k, 2))))...
                / (1-exp(1i*pi*(x_t-spatailf3(k, 1)))) / (1-exp(1i*pi*(y_t-spatailf3(k, 2))));
            else
            end
            spatailf3(k, 3) = sqrt(n_arr(3, 1)) * rp3(C3(k, 1), gr_n) / sqrt(rho) / lambdae / lambdamk;
        end