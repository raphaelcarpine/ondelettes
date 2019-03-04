for ModeNum=1:4
    
    QMode = Q(ModeNum) ;
    f_minMode = f_min(ModeNum) ;
    f_maxMode = f_max(ModeNum) ;
    nb_freq =       250 ;
    
    %% test input
    voies=1:4;
    X{1} = t;X{2}=t;X{3}=t;X{4}=t;
    Y{1} = y(:,1);Y{2} = y(:,3);Y{3} = y(:,5);Y{4} = y(:,7);
    %% ridge extract
    ridge=cell(length(voies),1);
    n_ridge=ridge;
    
    for C_voie=1:length(voies)
        ridge{C_voie} = RidgeExtract(X{C_voie},Y{voies(C_voie)},QMode,f_minMode,f_maxMode,nb_freq,...
            'MinModu',MinModu);
        n_ridge{C_voie}=length(ridge{C_voie}.val);
    end
    
    tf2 = 1.5*max([ridge{1}.time{1}(end),ridge{2}.time{1}(end),ridge{3}.time{1}(end),ridge{4}.time{1}(end)]);
    %% plot ridges freq
    figure('Name',sprintf('Test signal, mode %d',ModeNum),'NumberTitle','off')
    subplot(3,2,1)
    hold on
    ax=gca;
    for C_voie=1:length(voies)
        for C_r=1:n_ridge{C_voie}
            h=plot(ridge{C_voie}.time{C_r},ridge{C_voie}.freq{C_r},...
                'LineWidth',1.5,'DisplayName',sprintf('u_%d ridge %d',voies(C_voie),C_r));
            plot(ridge{C_voie}.time{C_r},ridge{C_voie}.freqraw{C_r},...
                'Color',get(h,'Color'),'LineStyle',':','HandleVisibility','off','LineWidth',1)
            ax.ColorOrderIndex=mod(ax.ColorOrderIndex-2,7)+1;
        end
    end
    
    xlabel('time [s]')
    ylabel('f(ridge) [Hz]')
    plot([0,tf2],[freq(ModeNum),freq(ModeNum)],'k--','DisplayName','exact')
    
    hold off
    axis([0,tf2,freq(ModeNum)*(1-PlotDiffExactInf),freq(ModeNum)*(1+PlotDiffExactSup)]);
    %% plot ridges abs
    subplot(3,2,2)
    hold on
    ax=gca;
    for C_voie=1:length(voies)
        for C_r=1:n_ridge{C_voie}
            h=plot(ridge{C_voie}.time{C_r},log(abs(ridge{C_voie}.val{C_r})),...
                'LineWidth',1.5,'DisplayName',sprintf('u_%d ridge %d',voies(C_voie),C_r));
            plot(ridge{C_voie}.time{C_r},log(abs(ridge{C_voie}.valraw{C_r})),...
                'Color',get(h,'Color'),'LineStyle',':','HandleVisibility','off','LineWidth',1);
            ax.ColorOrderIndex=mod(ax.ColorOrderIndex-2,7)+1;
        end
    end
    
    xlabel('time [s]')
    ylabel('log|CWT(ridge)|')
    hold off
    axis([0,tf2,-inf,+inf]);
    %% plot amor
    subplot(3,2,3)
    hold on
    ax=gca;
    for C_voie=1:length(voies)
        for C_r=1:n_ridge{C_voie}
            h=plot(ridge{C_voie}.time{C_r},ridge{C_voie}.inv2Q{C_r},...
                'LineWidth',1.5,'DisplayName',sprintf('u_%d ridge %d',voies(C_voie),C_r));
            plot(ridge{C_voie}.time{C_r},ridge{C_voie}.inv2Qraw{C_r},...
                'Color',get(h,'Color'),'LineStyle',':','HandleVisibility','off','LineWidth',1)
            ax.ColorOrderIndex=mod(ax.ColorOrderIndex-2,7)+1;
        end
    end
    
    xlabel('time [s]')
    ylabel('-diff log|CWT(ridge)|/(2\pif(ridge)) [1]')
    axis([0,tf2,xi(ModeNum)*(1-PlotDiffExactInf),xi(ModeNum)*(1+PlotDiffExactSup)]);
    
    plot([0,tf2],[xi(ModeNum),xi(ModeNum)],'k--','DisplayName','exact')
    hold off
    %% plot diff
    subplot(3,2,4)
    hold on
    ax=gca;
    for C_voie=1:length(voies)
        for C_r=1:n_ridge{C_voie}
            h=plot(ridge{C_voie}.time{C_r},ridge{C_voie}.bandwidth2{C_r},...
                'LineWidth',1.5,'DisplayName',sprintf('u_%d ridge %d',voies(C_voie),C_r));
            plot(ridge{C_voie}.time{C_r},ridge{C_voie}.bandwidthraw{C_r},...
                'Color',get(h,'Color'),'LineStyle',':','HandleVisibility','off','LineWidth',1)
            ax.ColorOrderIndex=mod(ax.ColorOrderIndex-2,7)+1;
        end
    end
    
    xlabel('time [s]')
    ylabel('-diff log|CWT(ridge)|/(2\pi) [Hz]')
    
    plot([0,tf2],freq(ModeNum)*[xi(ModeNum),xi(ModeNum)],'k--','DisplayName','exact')
    
    hold off
    
    axis([0,tf2,freq(ModeNum)*xi(ModeNum)*(1-PlotDiffExactInf),freq(ModeNum)*xi(ModeNum)*(1+PlotDiffExactSup)]);
    %%
    %% relative phase and amp
    ref0 = [4,1,1,2];
    ref = ref0(ModeNum);
    
    pha_ref = NaN(1,length(Y{ref}));
    [~,ind] = min(abs(X{ref}-ridge{ref}.time{1}(1)));
    pha_ref(ind:length(ridge{ref}.time{1})+ind-1) = ridge{ref}.pha{1};
    
    for C_voie=1:length(voies)
        for C_r=1:n_ridge{C_voie}
            ridge{C_voie}.relapha{C_r} = NaN(size(Y{C_voie}));
            ridge{C_voie}.amp{C_r} = NaN(size(Y{C_voie}));
            
            if X{ref} == X{C_voie}
                [~,ind] = min(abs(X{ref}-ridge{C_voie}.time{C_r}(1)));
                L_r = length(ridge{C_voie}.time{C_r});
                ridge{C_voie}.relapha{C_r}(ind:L_r+ind-1) = mod(ridge{C_voie}.pha{C_r} - pha_ref(ind:L_r+ind-1),2*pi);
            else
                error('pas le m?me ?chantillonnage')
            end
            
            ridge{C_voie}.amp{C_r}(ind:length(ridge{C_voie}.time{C_r})+ind-1) = abs(ridge{C_voie}.val{C_r});
            
        end
    end
    AmpRef = ridge{ref}.amp{1};
    for C_voie=1:length(voies)
        for C_r=1:n_ridge{C_voie}
            ridge{C_voie}.amp{C_r} = ridge{C_voie}.amp{C_r}./AmpRef;
        end
    end
    ax=polaraxes;
    Title = title('Amplitude-phase plot');
    hold on
    for C_voie=1:length(voies)
        for C_r=1:n_ridge{C_voie}
            h=polarplot(ridge{C_voie}.relapha{C_r},abs(ridge{C_voie}.amp{C_r}),...
                '-','LineWidth',1.5,'DisplayName',sprintf('u_%d ridge %d',voies(C_voie),C_r));
            I_1 = find(~isnan(abs(ridge{C_voie}.amp{C_r})),1);
            polarplot(ridge{C_voie}.relapha{C_r}(I_1),abs(ridge{C_voie}.amp{C_r}(I_1)),...
                'o','MarkerSize',6,'HandleVisibility','off','Color',get(h,'Color'),'MarkerFaceColor',get(h,'Color'));
            
            ax.ColorOrderIndex=mod(ax.ColorOrderIndex-2,7)+1;
        end
    end
    polarplot(angle(phi(:,ModeNum)),abs(phi(:,ModeNum)),...
        'ko','LineWidth',1.5,'HandleVisibility','off');
    polarplot(NaN,NaN,...
        'ko--','LineWidth',1,'DisplayName','exact');
    
    ax.OuterPosition = [.3249,-.0422,.3403,.46];
    Title.Position = [177,2.2,0];
    hold off
    Legend=legend('show');
    Legend.Position= [.6973,.0678,.1975,.2244];
    
    set(gcf,'Units','centimeters','Position',[0,0,14,16])
    saveas(gcf,fullfile(p,sprintf('SysLin4ddl_mode%d.png',ModeNum)))
    
    %% Mode Plot
    h=figure;
    h.Units = 'centimeters';
    h.Position = [0,0,[1,10]];
    hold on
    plot(0,0,'oblack','MarkerFaceColor','black')
    S(1)=plot(zeros(1,5),0:4,'-black');
    S(2)= plot(zeros(1,4),1:4,'ored','MarkerFaceColor','red');
    S(3)=plot(zeros(1,5),0:4,'oblack');
    
    axis([-1,1,-1,5])
    ax = gca;
    ax.Visible='off';
    
    ax.Position = ax.OuterPosition;
    
    
    for C_t = 0:1/50:1-1/50
        S(1).YData(2:5) = (1:4) + AmpFactor*real(transpose(phi(:,ModeNum))*exp(1i*2*pi*C_t));
        S(2).YData = (1:4) + AmpFactor*real(transpose(phi(:,ModeNum))*exp(1i*2*pi*C_t));
        S(3).YData(2:5) = (1:4) + AmpFactor*real(transpose(phi(:,ModeNum))*exp(1i*2*pi*C_t));
        
        drawnow
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        % Write to the GIF File
        if C_t == 0
            imwrite(imind,cm,fullfile(p,sprintf('SysLin4ddl_mode%d_ModeShape.gif',ModeNum)),'gif', 'Loopcount',inf,'DelayTime',1/40);
        else
            imwrite(imind,cm,fullfile(p,sprintf('SysLin4ddl_mode%d_ModeShape.gif',ModeNum)),'gif','WriteMode','append','DelayTime',1/40);
        end
    end
    close(h)
    
end