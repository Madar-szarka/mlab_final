%%% POSTROC
%%% ez a .m a hiba matrixon megy vegig a tesztelt ablakok szerint es rajzol
%%% idlefutas diagrammot az adott ablakokhoz tartozo legjobb a,b,c-vel 
%%% vagyis az a,b,c osszevissza valtozik
%%% (wss ablak darab szam) * (Tf2 ablak darab szam) abrat fog kesziteni
for i=1:win_WSS_ciklus_length
    win_WSS=win_WSS_ciklus_elements(i);
    win_Tf=win_WSS;
    
    T_WSS_smooth=ablak(T_WSS_A2,win_WSS);

    dT=T_WSS_smooth(2:end)-T_WSS_smooth(1:end-1);
    dt=time(2:end)-time(1:end-1);
    update_rate=mean(dt); % [s/consecutive sample]
    T_dot=dT./dt; 

    for j=1:win_Tf2_ciklus_length
        nnnnn=squeeze(error_matrix_RMS(:,:,:,i,j)); % 3D tombot emel ki adot i-j-re
        
        [vvvv,loc] = min(nnnnn(:));     % vvv=absolute min; loc=location (egy skalar index mert itt ki van teritve vektorra a matrix)
        [iiii,jjjj,kkkk] = ind2sub(size(nnnnn),loc); % legjobb illesztes 3 indexe a,b,c-re

        a_v = a_begin+(iiii-1)*a_lepes; %minimum ertekek, (NEM index)
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
        title (['winWSS = ',num2str(win_WSS*update_rate,formatspec_ablak),'   winTF2 = ',num2str(win_Tf2*update_rate,formatspec_ablak),...
            '   a_v_,_m_i_n = ',num2str(a_v,formatSpec),'   b_v_,_m_i_n = ',num2str(b_v,formatSpec),'   cRm_m_i_n = ',num2str(cRm_n1),...
                                 %{
%%% SZEDD KI HA MAR VAN MAE            '  MAE = ', num2str(MAE_matrix(iiii,jjjj,kkkk,i,j),'%.2f'),... %%%% MAE itt nem
MAE_min mivel a ciklus az RMS minimumat kersve fut
KI LEHETNE PROBALNI HOGY UGY FUSSON HOGY MAE SZERINT KERESI A LEGJOBBAT!!!!!!!!!!!!!!!!!!!!!
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

        saveas(gcf,['Temp_winWSS_',num2str(win_WSS*update_rate,'%.0f'),'_winTf2_',num2str(win_Tf2*update_rate,'%.0f'),...
            '_RMSE',num2str(error_matrix_RMS(iiii,jjjj,kkkk,i,j),'%.0f'),'.png']);

        close(fig);

    end
end


