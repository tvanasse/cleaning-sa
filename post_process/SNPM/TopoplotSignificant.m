function h = TopoplotSignificant(topodata,sigch,chanlocs,insidegoodch,titletext)
figure    
topoplot(topodata,chanlocs,'headrad',.56,'electrodes','on','maplimits','minmax','whitebk','on');
title(titletext);
colorbar
set(gca,'Xlim',[-.55 .55],'Ylim',[-.59 .59]);
electrodes.x=get(findobj(gca,'Marker','.'),'XData'); 
length(electrodes.x)
electrodes.y=get(findobj(gca,'Marker','.'),'YData'); 
electrodes.z=get(findobj(gca,'Marker','.'),'ZData');
delete(findobj(gca,'Marker','.'));

% if you use inside good, be careful to pull correct indicies 
if ~isequal(length(electrodes.x),length(insidegoodch))
    disp('error with plotting sig channels');
else
    if sum(sigch)>0
        [~,chi] = intersect(insidegoodch,sigch);
        hold on;h = scatter(electrodes.x(chi),electrodes.y(chi),electrodes.z(chi),...
            'filled','SizeData',80,'Cdata',[1 1 1],'MarkerEdgeColor',[0 0 0],'linewidth',.5);
    end
end
set(gcf,'InvertHardCopy','off','color','w');
