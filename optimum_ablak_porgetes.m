if exist('time','var')==0 %%% time mar a matlab-ban letrehozott ido ami az abszolut "t" meresbeli ido offset-el verzioja
time=t-(t(1));            %%% ha olyan kezelt merest toltunk be amiben tortent tordeles, akkor ott mar ezek aszerint letre lesznek hozva
v_ebc2_filt=filloutliers(MeanFASpeed_EBC2_CAN7,0,'mean');     

end

% win_WSS_ciklus_elements=5:5:50; % VT 10 30 50 szt komment 2
% win_WSS_ciklus_elements=[2 6 16];
%win_WSS_ciklus_elements=[10 30 60];
 win_WSS_ciklus_elements=10:10:100; 

win_WSS_ciklus_length=length(win_WSS_ciklus_elements);


%win_Tf2_ciklus_elements=[1:5:30,30:10:60];  % VT 10 60 90
%win_Tf2_ciklus_elements=[4 16 150]; % 05.14 ebed kozben
%win_Tf2_ciklus_elements=[10 30 60];
%win_Tf2_ciklus_elements=[0.5, 1, 2, 4 8 20 40];
win_Tf2_ciklus_elements=1:2:17;

win_Tf2_ciklus_length=length(win_Tf2_ciklus_elements);


            a_begin=0;            %%% a_begin=0;
            a_lepes=0.3/20*2;          %%% a_lepes=0.3;
            a_db=30;              %%% a_db=30;
            b_begin=0;            %%% b_begin=0;
            b_lepes=5/20;            %%% b_lepes=5;
            b_db=30;              %%% b_db=30;
            c_begin=200;          %%% c_begin=200;
            c_lepes=300;          %%% c_lepes=300;
            c_db=30;              %%% c_db=30;
            
            
            error_matrix=zeros(a_db,b_db,c_db,win_WSS_ciklus_length,win_Tf2_ciklus_length);
            error_matrix_RMS=zeros(a_db,b_db,c_db,win_WSS_ciklus_length,win_Tf2_ciklus_length);
            diff_matrix=zeros(a_db,b_db,c_db,win_WSS_ciklus_length,win_Tf2_ciklus_length);
            MAE_matrix=zeros(a_db,b_db,c_db,win_WSS_ciklus_length,win_Tf2_ciklus_length);
            
            min_error=2.000000000000000e+25;
            min_error_RMS=99999999999999;
            min_diff=2.000000000000000e+25;
            MAE=2.000000000000000e+25;
            
            
            aa=0;
            bb=0;
            cc=0;
            aa_min=0;
            bb_min=0;
            cc_min=0;
            aa_min_idx=0;
            bb_min_idx=0;
            cc_min_idx=0;
            win_WSS_min=0;
            win_Tf2_min=0;
            win_WSS_ciklus_idx=0;
            win_Tf2_ciklus_idx=0;
            win_WSS_ciklus_idx_min=0;
            win_Tf2_ciklus_idx_min=0;
            
            plot_rms_ind=1;
            plot_diff_ind=1;


for win_WSS=win_WSS_ciklus_elements  % win_Tf is figyelni
        win_WSS_ciklus_idx=win_WSS_ciklus_idx+1
        win_Tf2_ciklus_idx=0;
%win_WSS=20; % PAROS legyen fél ablak meret miatt% 20 sample szeles ablak
% % % % div_wss = win_WSS * ones(length(T_WSS_A2),1);
% % % % div_wss(1:win_WSS) = (1:win_WSS); % elso ablaknyi elemnel nem ablakmerettel osztunk
% % % % sum_wss = zeros(length(T_WSS_A2),1);

T_WSS_smooth=ablak(T_WSS_A2,win_WSS);


    for win_Tf2_factor= win_Tf2_ciklus_elements      % win_Tf2 / win_WSS ratio | korabban 30*
        win_Tf2_ciklus_idx=win_Tf2_ciklus_idx+1
        aa=0;
        win_Tf=win_WSS;
        win_Tf2=win_Tf2_factor*win_WSS;

%             cRm_n1=5505;        %%%  cRm_n1=5505;
%             a_v=7/80/20;           %%%  a_v=7/80;
%             b_v=1.25/20;           %%%  b_v=1.25;
            
            
            % Tf_n1_smooth ablak merete ha NEM win_Tf=win_WSS
            % gradT_lim meg megfontolando
            
            
% %             oszto =  win_Tf * ones(length(T_WSS_smooth)-1,1); % azert egyel rovidebb mert majd derivalt mennyiseg simitasahoz fogjuk felhasznalni
% %             oszto(1: win_Tf) = (1: win_Tf);
% %             osszeg = zeros(length(T_WSS_smooth)-1,1);
% %             
% %             oszto2 =  win_Tf2 * ones(length(T_WSS_smooth)-1,1); % azert egyel rovidebb mert majd derivalt mennyiseg simitasahoz fogjuk felhasznalni
% %             oszto2(1: win_Tf2) = (1: win_Tf2);
%           
% %             % T_WSS gradiens
            dT=T_WSS_smooth(2:end)-T_WSS_smooth(1:end-1);
            dt=time(2:end)-time(1:end-1);
            T_dot=dT./dt; 
            update_rate=mean(dt); % [s/consecutive sample]

          
            
            
            
  %%% kiszedve amig nem porgetem 
            for a_v=a_begin:a_lepes:(a_lepes*(a_db-1)+a_begin)
                aa=aa+1
                bb=0;
                for b_v=b_begin:b_lepes:(b_lepes*(b_db-1)+b_begin)
                    bb=bb+1;
                    cc=0;
                    
                    fileID = fopen('result_0515_ejjel_win_WSS_150_50kent_500_win_Tffactor_0dot5_x2_8ig.txt','at');
                    fprintf(fileID,'a_v index=%3d | b_v index=%3d | win_WSS_ciklus_idx=%3d | win_Tf2_ciklus_idx=%3d || minimum error=%12.0f | RMS error=%03d | a_v min/index=%11.5f / %3d | b_v min/index=%11.5f / %3d | CRm_n1 min/index=%11.5f / %3d | win_WSS_min/index=%4d / %4d | win_Tf2_min/index=%4d / %4d    ---> %s \n',...
                                    aa,          bb,                win_WSS_ciklus_idx,      win_Tf2_ciklus_idx,    min_error,               min_error_RMS,  aa_min,aa_min_idx,       bb_min,bb_min_idx,        cc_min,cc_min_idx,              win_WSS_min, win_WSS_ciklus_idx_min , win_Tf2_min,win_Tf2_ciklus_idx_min, datetime);
                    fclose(fileID);
                    
                    for cRm_n1=c_begin:c_lepes:(c_lepes*(c_db-1)+c_begin) % c=490 j/kg*K | m:=10kg | R=x/A_k ; 
                        cc=cc+1;
                        
                        R_cs_n1=a_v*v_ebc2_filt+b_v;
                        Tf_n1 = (1+R_cs_n1(2:end)).*T_WSS_smooth(2:end) + cRm_n1*T_dot(:) - R_cs_n1(2:end).*T_Ambient_CANSAS(2:end);

                        % eredmeny simitasa ablakkal
                        Tf_n1_smooth = ablak(Tf_n1,win_Tf); 
                        
                        % nagy ugrasok kiszedese gradiens alapjan
                        gradT_lim= 15; % Celsius/sec
                        dt_Tf=time(2:end)-time(1:end-1);
                        dTf_n1_smooth=Tf_n1_smooth(2:end)-Tf_n1_smooth(1:end-1); % ez mar ketto elemmel rovidebb mint az eredeti mert signalok


                        Tf_n1_smooth_gradfilt=Tf_n1_smooth;
                       % kiugro=4*ones(length(Tf_n1_smooth),1); %4 ha nem nyult hozza, 1 ha
                        for i=1:length(Tf_n1_smooth)-1
                            if  (Tf_n1_smooth_gradfilt(i+1,1) - Tf_n1_smooth_gradfilt(i,1)) / dt_Tf(i) > gradT_lim % ha extrem emelkedik
                                    Tf_n1_smooth_gradfilt(i+1,1) = Tf_n1_smooth_gradfilt(i,1) + gradT_lim * dt_Tf(i);
                                   % kiugro(i,1)=1;
                            elseif (Tf_n1_smooth_gradfilt(i+1,1) - Tf_n1_smooth_gradfilt(i,1)) / dt(i) < -gradT_lim % ha extrem csokken
                                    Tf_n1_smooth_gradfilt(i+1,1) = Tf_n1_smooth_gradfilt(i,1) - gradT_lim * dt_Tf(i);
                                  %  kiugro(i,1)=1; 
                            else
                              %  kiugro(i,1)=0;
                            end
                        end
                        
                        
                        Tf_n1_smoother = ablak(Tf_n1_smooth_gradfilt,win_Tf2); 
                        
                        
                        error=(Tf_n1_smoother-T_Brake_disc_A2(2:end)).^2; 
                        sum_error=sum(error); 
                        error_matrix(aa,bb,cc,win_WSS_ciklus_idx,win_Tf2_ciklus_idx)=sum_error;
                        error_RMS=sqrt(sum_error/(length(time)-1)); 
                        error_matrix_RMS(aa,bb,cc,win_WSS_ciklus_idx,win_Tf2_ciklus_idx)=error_RMS;
                        
                        difff=max(abs(Tf_n1_smoother-T_Brake_disc_A2(2:end)));
                        diff_matrix(aa,bb,cc,win_WSS_ciklus_idx,win_Tf2_ciklus_idx)=difff;
                        MAE=mean(abs(Tf_n1_smoother-T_Brake_disc_A2(2:end))); %%% Mean Absolute Error
                        MAE_matrix(aa,bb,cc,win_WSS_ciklus_idx,win_Tf2_ciklus_idx)=MAE;
                        
                        if sum(error) < min_error
                            min_error=sum(error);
                            min_error_RMS=error_RMS;
                            aa_min_idx=aa;
                            bb_min_idx=bb;
                            cc_min_idx=cc;
                            aa_min=a_v;
                            bb_min=b_v;
                            cc_min=cRm_n1;
                            win_WSS_ciklus_idx_min=win_WSS_ciklus_idx;
                            win_Tf2_ciklus_idx_min=win_Tf2_ciklus_idx;
                            win_WSS_min=win_WSS;
                            win_Tf2_min=win_Tf2;
                            
                            %%%% ez a plot egyre jobban illeszkedo gorbet
                            %%%% rajzol, ahogy fut az 5D ciklus, de
                            %%%% következetesseg nincs benne mert
                            %%%% tetszeloges index ugrasok lehetnek a
                            %%%% talalt legjobb illesztesek kozott
                            
                            fig=figure('units','normalized','Position',[0.1 0.05 0.8 0.85],'Visible','off'); % visible off hogy ne ugorjon fel az ablak
                            plot(time,T_WSS_A2,'k-',time,T_Brake_disc_A2,'b-',time(2:end,1),Tf_n1_smoother,'g-','LineWidth',1);
                            lgd=legend("Meas. sensor", "Meas. brake", "Calc. brake ",'Location','northwest');
                            lgd.FontSize = 16.0;

                            ylim([0 800]);

                            formatSpec = '%.3f';
                            title (['winWSS = ',num2str(win_WSS*update_rate),'   winTF2 = ',num2str(win_Tf2*update_rate),'   a_v = ',num2str(a_v,formatSpec),'   b_v = ',num2str(b_v,formatSpec),'   cRm = ',num2str(cRm_n1),'   max(abs(Diff)) = ',num2str(difff),'   RMSE = ',num2str(min_error_RMS)],'FontSize',15);

                            
                            ax = gca;
                            ax.GridAlpha = 0.3;
                            grid on;
                            xAX = get(gca,'XAxis');
                            set(xAX,'FontSize', 20);
                            yAX = get(gca,'YAxis');
                            set(yAX,'FontSize', 20);
                            xlabel('Time [s]', 'FontSize',20);
                            ylabel('Temperature[°C]','FontSize',20);

                            saveas(gcf,['T_vs_t_RMS_win_WSS',num2str(win_WSS*update_rate,'%.0f'),'_winTf2_',num2str(win_Tf2*update_rate,'%.0f'),'__a_v',num2str(a_v,'%.0f'),'_b_v',num2str(b_v,'%.0f'),'_cRm',num2str(cRm_n1,'%.0f'),'_',num2str(plot_rms_ind),'.png']);
                            frame = getframe(fig);
                            im{plot_rms_ind} = frame2im(frame);
                            close(fig);
  

                            
                            plot_rms_ind=plot_rms_ind+1;
                          end
                          
                         if difff < min_diff
                            min_diff=difff;


                            fig=figure('units','normalized','Position',[0.1 0.05 0.8 0.85],'Visible','off');
                            plot(time,T_WSS_A2,'k-',time,T_Brake_disc_A2,'b-','LineWidth',1);
                            hold on;
                            plot(time(2:end,1),Tf_n1_smoother,'g-','LineWidth',1);
                            lgd=legend("Meas. sensor", "Meas. brake","Calc. brake ",'Location','northwest');
                            lgd.FontSize = 16.0;

                            ylim([0 1000]);

                            formatSpec = '%.3f';
                            title (['winWSS = ',num2str(win_WSS*update_rate),'   winTF2 = ',num2str(win_Tf2*update_rate),'   a_v = ',num2str(a_v,formatSpec),'   b_v = ',num2str(b_v,formatSpec),'   c_v = ',num2str(cRm_n1),'   max(abs(Diff)) = ',num2str(difff),'   RMS = ',num2str(min_error_RMS)],'FontSize',15);


                            ax = gca;
                            ax.GridAlpha = 0.3;
                            grid on;
                            xAX = get(gca,'XAxis');
                            set(xAX,'FontSize', 20);
                            yAX = get(gca,'YAxis');
                            set(yAX,'FontSize', 20);
                            xlabel('Time [s]', 'FontSize',20);
                            ylabel('Temperature[°C]','FontSize',20);

                            saveas(gcf,['T_vs_t_diff_win_WSS',num2str(win_WSS*update_rate,'%.0f'),'_winTf2_',num2str(win_Tf2*update_rate,'%.0f'),'__a_v',num2str(a_v,'%.0f'),'_b_v',num2str(b_v,'%.0f'),'_cRm',num2str(cRm_n1,'%.0f'),'_',num2str(plot_diff_ind),'.png']);
                            frame = getframe(fig);
                            im_diff{plot_diff_ind} = frame2im(frame);
                            close(fig);

                            plot_diff_ind=plot_diff_ind+1;

                        end
                          
                       end
                    end
                end
    end
   %filename_error_mtx_tolig=sprintf('%d_%d',win_WSS(1),win_WSS_ciklus_elements(win_WSS_ciklus_idx));   
   filenev=sprintf('error_matrix_%02d_%02d_%02d_win_WSS_%03d_%03d',day(datetime),hour(datetime),minute(datetime),win_WSS_ciklus_elements(1),win_WSS_ciklus_elements(win_WSS_ciklus_idx));
   save(filenev,'error_matrix');
   
   filenev_RMS=sprintf('RMS_error_matrix_%02d_%02d_%02d_win_WSS_%03d_%03d',day(datetime),hour(datetime),minute(datetime),win_WSS_ciklus_elements(1),win_WSS_ciklus_elements(win_WSS_ciklus_idx));
   %%%filenev_RMS=sprintf('RMS_error_matrix_%02d_%02d_%02d_win_WSS_%03d',day(datetime),hour(datetime),minute(datetime),win_WSS_ciklus_elements(win_WSS_ciklus_idx));
   save(filenev_RMS,'error_matrix_RMS'); 
   
   filenev_MAE=sprintf('MAE_matrix_%02d_%02d_%02d_win_WSS_%03d_%03d',day(datetime),hour(datetime),minute(datetime),win_WSS_ciklus_elements(1),win_WSS_ciklus_elements(win_WSS_ciklus_idx));
   save(filenev_MAE,'MAE_matrix'); 
   
end


   filenev_workspace=sprintf('Workspace_finish_%02d_%02d_%02d_win_WSS_%03d_%03d',day(datetime),hour(datetime),minute(datetime),win_WSS_ciklus_elements(1),win_WSS_ciklus_elements(win_WSS_ciklus_idx));
   save(filenev_workspace);
% % %  %%mukodo
% % % filename_gif = "testoptimum_balak_porgetes_0515.gif"; % Specify the output file name
% % % for i_gif = 1:plot_rms_ind-1
% % %     [A,map] = rgb2ind(im{i_gif},256);
% % %     if i_gif == 1
% % %         imwrite(A,map,filename_gif,"gif","LoopCount",Inf,"DelayTime",1);
% % %     else
% % %         imwrite(A,map,filename_gif,"gif","WriteMode","append","DelayTime",1);
% % %     end
% % % end



