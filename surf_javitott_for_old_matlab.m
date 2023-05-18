
error_matrix_RMS_3D=squeeze(error_matrix_RMS(:,:,:,win_WSS_ciklus_idx_min,win_Tf2_ciklus_idx_min)); % a kb legjobb illeszteshez tartozo ablakmeret indexet beirni

% ellenorizni hogy mik a szelsoerteke az RMSE-nek es az alapjan loni be az
% abra hatarokat
max(max(max(error_matrix_RMS_3D)))
min(min(min(error_matrix_RMS_3D)))



for i=1:a_db
    fig=figure('units','normalized','Position',[0.1 0.0 0.6 0.93],'Visible','off');
    nnnnn=squeeze(error_matrix_RMS_3D(i,:,:)); % i.-edik lapot (2D metszet) kiemeli | !! b es c sorrendje ahogy kiemeli fontos
    xxxxx = [c_begin (c_begin+c_lepes * (c_db-1))]; % elso es utolso elem (2 db)
    yyyyy = [b_begin (b_begin+b_lepes * (b_db-1))]; % elso es utolso elem (2 db)

    imagesc(xxxxx,yyyyy,nnnnn)
    colorbar
    cbar = colorbar;
    cbar.Label.String = 'RMSE';
%     cbar.TickLabels={'cold','warm','hot'};
    cbar.YTick = [50 300 660 ]; % ha     ax.CLim=[** ***]; modosul akkor utana kellhet huzni
    cbar.YTickLabel = {'50', '300', '660'}; % szigoruan koveti cbar.YTick  ertekeket
    
    ax = gca;
    
    ax.ColorScale = 'log';
    ax.CLim=[50 660];
    
    axis square

    xlabel('cRm', 'FontSize',18);
    ylabel('b_v','FontSize',18);
    ax.FontSize = 15;
    %---------------------------------------------

    [AAAmin,loc] = min(nnnnn(:));
    [AAAi,AAAj] = ind2sub(size(nnnnn),loc);


    hold on
    plot((c_begin+(AAAj-1)*c_lepes),(b_begin+(AAAi-1)*b_lepes), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
    formatSpec = '%.3f';
    title (['winWSS_m_i_n = ',num2str(win_WSS_min),'   winTf2_m_i_n = ',num2str(win_Tf2_min),'   a_v = ',num2str(a_begin+(i-1)*a_lepes,formatSpec),'   b_v_,_m_i_n = ',num2str(b_begin+(AAAi-1)*b_lepes,formatSpec),'   cRm_m_i_n = ',num2str(c_begin+(AAAj-1)*c_lepes),'   RMSE_m_i_n = ',num2str(AAAmin,'%.2f')],'FontSize',15);
    %subtitle(['winWSS_m_i_n = ',num2str(win_WSS_min),'   winTf2_m_i_n = ',num2str(win_Tf2_min)],'FontSize',10);
    
    ax.YDir = 'normal';


    saveas(gcf,['a_surf_',num2str(i),'.png']);

    %exportgraphics(gcf,'a_testAnimated.gif','Append',true);

    close(fig);
    % pause(1) gif create miatt
end



for i=1:b_db
    fig=figure('units','normalized','Position',[0.1 0.0 0.6 0.93],'Visible','off');
    nnnnn=squeeze(error_matrix_RMS_3D(:,i,:));
    xxxxx = [c_begin (c_begin+c_lepes * (c_db-1))];
    yyyyy = [a_begin (a_begin+a_lepes * (a_db-1))];

     imagesc(xxxxx,yyyyy,nnnnn)
    colorbar
    cbar = colorbar;
    cbar.Label.String = 'RMSE';
    cbar.YTick = [50 300 660 ]; % ha     ax.CLim=[** ***]; modosul akkor utana kellhet huzni
    cbar.YTickLabel = {'50', '300', '660'}; % szigoruan koveti cbar.YTick  ertekeket
    
    
    ax = gca;
    
    ax.ColorScale = 'log';
    ax.CLim=[50 660];
    
    axis square

    xlabel('cRm', 'FontSize',18);
    ylabel('a_v','FontSize',18);
    ax.FontSize = 15;

    %---------------------------------------------

    [AAAmin,loc] = min(nnnnn(:));
    [AAAi,AAAj] = ind2sub(size(nnnnn),loc);


    hold on
    plot((c_begin+(AAAj-1)*c_lepes),(a_begin+(AAAi-1)*a_lepes), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
    formatSpec = '%.3f';
    title (['b_v = ',num2str(b_begin+(i-1)*b_lepes,formatSpec),'   a_v_,_m_i_n = ',num2str(a_begin+(AAAi-1)*a_lepes,formatSpec),'   cRm,_m_i_n = ',num2str(c_begin+(AAAj-1)*c_lepes),'   RMSE_m_i_n = ',num2str(AAAmin,'%.2f')],'FontSize',15);

    ax.YDir = 'normal';
    

    saveas(gcf,['b_surf_',num2str(i),'.png']);

  %  exportgraphics(gcf,'b_testAnimated.gif','Append',true);

    close(fig);
  %  pause(1)
end


for i=1:c_db
    fig=figure('units','normalized','Position',[0.1 0.0 0.6 0.93],'Visible','off');
    nnnnn=squeeze(error_matrix_RMS_3D(:,:,i));
    xxxxx = [b_begin (b_begin+b_lepes * (b_db-1))];
    yyyyy = [a_begin (a_begin+a_lepes * (a_db-1))];

    imagesc(xxxxx,yyyyy,nnnnn)
    colorbar
    cbar = colorbar;
    cbar.Label.String = 'RMSE';
    cbar.YTick = [50 300 660 ]; % ha     ax.CLim=[** ***]; modosul akkor utana kellhet huzni
    cbar.YTickLabel = {'50', '300', '660'}; % szigoruan koveti cbar.YTick  ertekeket
    
    ax = gca;
    
    ax.ColorScale = 'log';
    ax.CLim=[50 660];
    
    axis square

    xlabel('b_v', 'FontSize',18);
    ylabel('a_v','FontSize',18);
    ax.FontSize = 15;

    %---------------------------------------------

    [AAAmin,loc] = min(nnnnn(:));
    [AAAi,AAAj] = ind2sub(size(nnnnn),loc);


    hold on
    plot((b_begin+(AAAj-1)*b_lepes),(a_begin+(AAAi-1)*a_lepes), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
    formatSpec = '%.3f';
    title (['cRm = ',num2str(c_begin+(i-1)*c_lepes,formatSpec),'   a_v_,_m_i_n = ',num2str(a_begin+(AAAi-1)*a_lepes,formatSpec),'   b_v_,_m_i_n = ',num2str(b_begin+(AAAj-1)*b_lepes),'   RMSE_m_i_n = ',num2str(AAAmin,'%.2f')],'FontSize',15);

    ax.YDir = 'normal';
    

    saveas(gcf,['c_surf_',num2str(i),'.png']);

   % exportgraphics(gcf,'c_testAnimated.gif','Append',true);

    close(fig);
%    pause(1)
end
