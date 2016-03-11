name_set=textread('grk2011_Header.csv','%s');
period=1;
mi=data.InormByDist(:,:,1,period);
th(:,:)=data.SigThreshTvsIzero(:,:,period);
Tz_result(:,:,:)=data.TvsIzero(:,:,:,period);
Tz_result(:,:,1)=mi./mi(:,:);
x=(0:0.5:18);
[a,~]=size(name_set);


for j=2
    fig=figure(j);
    data_var=name_set{j,1};
    for i=10
        y_ij(:,1)=Tz_result(j,i,:);
        y_ji(:,1)=Tz_result(i,j,:);
        th_ij(1,1)=data.SigThreshInormByDist(i,j,period)./mi(i,j);
        th_ji(1,1)=data.SigThreshInormByDist(j,i,period)./mi(j,i);
        for k=2:37
            th_ij(k,1)=th(j,i);
            th_ji(k,1)=th(i,j);
        end        
        x_g=strcat(data_var,'>',name_set{i,1});
        y_g=strcat(name_set{i,1},'>',data_var);
        %subplot(4,4,i);        
        plot(x,y_ij,'-k','LineWidth',2)
        hold
        plot(x,y_ji,'-','color',[0.5 0.5 0.5],'LineWidth',2);
        plot(x,th_ij,'-.k','LineWidth',2);
        plot(x,th_ji,'-.','color',[0.5 0.5 0.5],'LineWidth',2);
        xlim([0 18]);
        ylabel('Tz=T/I');
        xlabel('Time lag[h]')
        hleg=legend(x_g,y_g);
        set(hleg,'FontSize',7)
        set(hleg,'Color','none')
        style=hgexport('factorystyle');
        style.Bounds='Tight';
        hgexport(gcf,'-Clipboard',style,'ApplyStyle',true);
        saveas(fig,strcat('canonical_type','_',data_var,'_',name_set{i,1},'.png'));
        hold 
    end
end

        
        
        
