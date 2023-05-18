%%% POSTROC
%%% szig monoton javuloan megy vegig az 5D for ciklus porgetese
%%% sorrendjeben

current_error=10e25;
plot_rms_ind_post=0;

for i=1:win_WSS_ciklus_length
    for j=1:win_Tf2_ciklus_length
        for iiii=1:a_db
           for jjjj=1:b_db
              for kkkk=1:c_db
              
                  if error_matrix_RMS(iiii,jjjj,kkkk,i,j) < current_error
                    
                    plot_rms_ind_post=plot_rms_ind_post+1;  
                    current_error=error_matrix_RMS(iiii,jjjj,kkkk,i,j);

                    win_WSS=win_WSS_ciklus_elements(i);
                    win_Tf=win_WSS;

                    T_WSS_smooth=ablak(T_WSS_A2,win_WSS);

                    dT=T_WSS_smooth(2:end)-T_WSS_smooth(1:end-1);
                    dt=time(2:end)-time(1:end-1);
                    update_rate=mean(dt); % [s/consecutive sample]
                    T_dot=dT./dt; 

                    a_v = a_begin+(iiii-1)*a_lepes; 
                    b_v = b_begin+(jjjj-1)*b_lepes;
                    cRm_n1 = c_begin+(kkkk-1)*c_lepes;


                    win_Tf2_factor= win_Tf2_ciklus_elements(j);
                    win_Tf2=win_Tf2_factor*win_WSS;

                    R_cs_n1=a_v*v_ebc2_filt+b_v;
                    Tf_n1 = (1+R_cs_n1(2:end)).*T_WSS_smooth(2:end) + cRm_n1*T_dot(:) - R_cs_n1(2:end).*T_Ambient_CANSAS(2:end); %  

                    % eredmeny simitasa ablakkal
                    Tf_n1_smooth = ablak(Tf_n1,win_Tf);

                    % nagy ugrasok kiszedese gradiens alapjan
                    gradT_lim= 15; % Celsius/sec
                    dt_Tf=time(2:end)-time(1:end-1);
                    dTf_n1_smooth=Tf_n1_smooth(2:end)-Tf_n1_smooth(1:end-1); % ez mar ketto elemmel rovidebb mint az eredeti mert signalok


                    Tf_n1_smooth_gradfilt=Tf_n1_smooth;
                    for ii=1:length(Tf_n1_smooth)-1
                        if  (Tf_n1_smooth_gradfilt(ii+1,1) - Tf_n1_smooth_gradfilt(ii,1)) / dt_Tf(ii) > gradT_lim % ha extrem emelkedik
                            Tf_n1_smooth_gradfilt(ii+1,1) = Tf_n1_smooth_gradfilt(ii,1) + gradT_lim * dt_Tf(ii);
                        elseif (Tf_n1_smooth_gradfilt(ii+1,1) - Tf_n1_smooth_gradfilt(ii,1)) / dt(ii) < -gradT_lim % ha extrem csokken
                            Tf_n1_smooth_gradfilt(ii+1,1) = Tf_n1_smooth_gradfilt(ii,1) - gradT_lim * dt_Tf(ii);
                        end
                    end


                    Tf_n1_smoother = ablak(Tf_n1_smooth_gradfilt,win_Tf2);


                    fig=figure('units','normalized','Position',[0.1 0.05 0.8 0.85],'Visible','off');
                    plot(time,T_WSS_A2,'k-',time,T_Brake_disc_A2,'b-',time(2:end,1),Tf_n1_smoother,'g-','LineWidth',1);
                    lgd=legend("Meas. sensor", "Meas. brake", "Calc. brake ",'Location','northwest');
                    lgd.FontSize = 16.0;

                    ylim([0 600]);

                    formatSpec = '%.3f';
                    formatspec_ablak='%.1f';
                    title (['winWSS = ',num2str(win_WSS*update_rate,formatspec_ablak),...
                        '   winTF2 = ',num2str(win_Tf2*update_rate,formatspec_ablak),...
                        '   a_v = ',num2str(a_v,formatSpec),...
                        '   b_v = ',num2str(b_v,formatSpec),...
                        '   cRm = ',num2str(cRm_n1),...
                        '   max(abs(Diff)) = ',num2str(diff_matrix(iiii,jjjj,kkkk,i,j),'%.2f'),...
                     %{
%%% SZEDD KI HA MAR VAN MAE            '  MAE = ', num2str(MAE_matrix(iiii,jjjj,kkkk,i,j),'%.2f'),...
                       %} 
                        '   RMSE_m_i_n = ',num2str(error_matrix_RMS(iiii,jjjj,kkkk,i,j),'%.2f')],'FontSize',15);

                    ax = gca;
                    ax.GridAlpha = 0.3;
                    grid on;
                    xAX = get(gca,'XAxis');
                    set(xAX,'FontSize', 20);
                    yAX = get(gca,'YAxis');
                    set(yAX,'FontSize', 20);
                    xlabel('Time [s]', 'FontSize',20);
                    ylabel('Temperature [°C]','FontSize',20);
                    
                    frame = getframe(fig);
                    im_post{plot_rms_ind_post} = frame2im(frame);

                    saveas(gcf,['T_vs_t_for5D_RMS_winWSS_',num2str(win_WSS*update_rate,'%.0f'),...
                        '_winTf2_',num2str(win_Tf2*update_rate,'%.0f'),'_a',num2str(a_v,'%.1f'),...
                        '_b',num2str(b_v,'%.1f'),'_c',num2str(cRm_n1,'%4d'),'_RMSE',num2str(error_matrix_RMS(iiii,jjjj,kkkk,i,j),'%.0f'),'.png']);
                    close(fig);
                  end
              end
           end
        end
    end
end

%%%%% GIF gyartas 1.
filename_gif_post = (['T_vs_t_for5D_RMS_winWSSmax_',num2str(win_WSS*update_rate,'%.0f'),...
                        '_winTf2max_',num2str(win_Tf2*update_rate,'%.0f'),'_amax',num2str(a_v,'%.1f'),...
                        '_bmax',num2str(b_v,'%.1f'),'_cmax',num2str(cRm_n1,'%4d'),'_RMSE',num2str(current_error,'%.0f'),'.gif']); % Specify the output file name
%filename_gif_post = (['testoptimum_ablak_porgetes_post_',num2str(day(datetime)),num2str(hour(datetime)),num2str(minute(datetime)),'.gif']); % Specify the output file name
for i_gif_post = 1:plot_rms_ind_post-1
    [A,map] = rgb2ind(im_post{i_gif_post},256);
    if i_gif_post == 1
        imwrite(A,map,filename_gif_post,"gif","LoopCount",Inf,"DelayTime",0.1);
    else
        imwrite(A,map,filename_gif_post,"gif","WriteMode","append","DelayTime",0.1);
    end
end


%%%%%% GIF gyartas 2.
%%%% result workspace rajzokbol
% %filename_gif = (['testoptimum_ablak_porgetes_',num2str(day(datetime)),num2str(hour(datetime)),num2str(minute(datetime)),'.gif']); % Specify the output file name
% % for i_gif = 1:plot_rms_ind-1
% %     [A,map] = rgb2ind(im{i_gif},256);
% %     if i_gif == 1
% %         imwrite(A,map,filename_gif,"gif","LoopCount",Inf,"DelayTime",0.1);
% %     else
% %         imwrite(A,map,filename_gif,"gif","WriteMode","append","DelayTime",0.1);
% %     end
% % end
