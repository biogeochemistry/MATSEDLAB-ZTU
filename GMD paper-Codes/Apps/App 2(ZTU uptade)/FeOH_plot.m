function FeOH_plot
    close all;
    global SimValues;
    x = linspace(0,15,511);             %defining spatial domain, must be the same as for sol
    conc = zeros(4,301);                  %initializing the matrix for partial derivative
    
    for xout = [0,1,2,3]
        for t1 = 1:101
            [u,~] = pdeval(0,x,SimValues{2,1},xout*5);  %calculating u at each t1
            conc(xout+1,t1) = u;
        end

        s_index = 2;

        for t2 = 102:301
            [u,~] = pdeval(0,x,SimValues{2,s_index},xout);  %calculating dudx at each t2
            conc(xout+1,t2) = u;
            s_index = s_index+1;
        end
    
        subplot(4,2,xout*2+1);
        plot([0:99,100:0.5:200],conc(xout+1,:),'LineWidth',1.5);
        hold on;
        set(gcf,'Position',get(0,'ScreenSize'));    %maximize the figure
        set(gcf,'PaperPosition',[3,5,24,11]);       %set the dimensions on paper when the plots are printed
        set(gca,'FontSize',14);                     %set the font size

        title('FeOH', 'FontSize',16,'FontWeight','bold');
        xlabel('Time (yr)','FontSize',12);
        ylabel([num2str(xout*5),' cm (umol/g)'],'FontSize',12);

        subplot(4,2,xout*2+2);
        plot(120:0.5:200,conc(xout+1,141:301),'LineWidth',1.5);
        hold on;
        set(gcf,'Position',get(0,'ScreenSize'));    %maximize the figure
        set(gcf,'PaperPosition',[3,5,24,11]);       %set the dimensions on paper when the plots are printed
        set(gca,'FontSize',14);                     %set the font size
        
        title('FeOH', 'FontSize',16,'FontWeight','bold');
        xlabel('Time (yr)','FontSize',12);
        ylabel([num2str(xout*5),' cm (umol/g)'],'FontSize',12);
    end
end