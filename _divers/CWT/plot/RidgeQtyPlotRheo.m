function [AxOut] = RidgeQtyPlotRheo(ridge,QtyX,QtyY,NameX,NameY,ScaleX,ScaleY,LimX,LimY)
%% Quantity : time, val, freq, diff, amor, freq2, pha

hold on
AxOut=gca;
NbPlot=1;
for C_voie=1:length(ridge)
    n_ridge=length(ridge{C_voie}.val);         %AxOut.ColorOrderIndex
    
    for C_r=1:n_ridge
        x       = eval(['ridge{C_voie}.',QtyX,'{C_r}']);
        y       = eval(['ridge{C_voie}.',QtyY,'{C_r}']);
        ybrut   = eval(['ridge{C_voie}.',QtyY,'raw','{C_r}']);
        
        if any(imag(x)~=0)
            x=abs(x);
        elseif any(imag(y)~=0)
            y=abs(y);
            ybrut=abs(ybrut);
        end
        if strcmp(ScaleX,'log')
            x=log(x);
        end
        if strcmp(ScaleY,'log')
            y=log(y);
            ybrut=log(ybrut);
        end
        
        h= plot(x,ybrut,...
            'LineStyle',':','HandleVisibility','off','LineWidth',1);
        plot(x,y,...
            'Color',get(h,'Color'),'LineWidth',1,'HandleVisibility','off');
        if NbPlot<8%8
            MarkerStr = 'o';
        elseif NbPlot<15%15
            MarkerStr = 'p';
        else
            MarkerStr = '*';
        end
        I1 = find(~isnan(y),1);
        I2 = find(~isnan(y),1,'last');
        plot([x(I1),x(I2)],[y(I1),y(I2)],...
            MarkerStr,'Color',get(h,'Color'),'MarkerSize',6,'HandleVisibility','off');
        
        
        
        %if C_voie==length(voies)||(voies(C_voie+1)~=voies(C_voie))
        plot(NaN,NaN,strcat(MarkerStr,'-'),'LineWidth',1,'Color',get(h,'Color'),...
            'DisplayName',sprintf('ridge %d',C_r))
        
        % end
        % AxOut.ColorOrderIndex
        %AxOut.ColorOrderIndex=mod(NbPlot,7)+1;
        NbPlot=NbPlot+1; %% ATTENTION SI PLUSIEURS RIDGES SUR MEME VOIE
        % end
        
        AxOut.ColorOrderIndex=mod(NbPlot-1,7)+1;
        
    end
end

%legend('show','location','best')
xlabel(NameX)
ylabel(NameY,'Interpreter','latex','FontSize',13)
% AxOut.XScale = ScaleX;
% AxOut.YScale = ScaleY;
AxOut.YGrid = 'on';
AxOut.XLim = LimX;
AxOut.YLim = LimY;
hold off
end