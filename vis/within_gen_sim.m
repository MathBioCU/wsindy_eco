%%%% view sims performed during algorithm to estimat Y(T)

% legend({'data','learned'},'location','westoutside')
legs = {'$S_{','$P_{1,','$P_{2,','$\nu_{1,','$\nu_{2,'};
subsamp = 1;
X_sub = find(diff(train_inds)==1);
ind=1;
for j=[1 2]
    for i=1:5
        if ismember(X_sub(i),X_in)
            figure(30+ind);
            set(gcf,'pos',[1254         589         308         299])
            plot(Y_ns{i,2}(1:subsamp:end),Y_ns{i,1}(1:subsamp:end,j),'b-',...
                train_time{i},nY(j)*Y_train{X_in==X_sub(i)}(:,j),'r.',...
                train_time{i}(end)*[1;1],[min(Y_ns{i,1}(:,j));max(Y_ns{i,1}(:,j))],'m--',...
                'linewidth',3,'markersize',14)
            xlabel('$t$','interpreter','latex')
            legend([legs{j},num2str(i),'}$'],'interpreter','latex')
            set(gca,'fontsize',24,'ticklabelin','latex',...
                'Yscale','linear','Xtick',0:ceil(yearlength/4):yearlength)
            grid on
            xlim([0 Y_ns{i,2}(end)])
            ylim([min([Y_ns{i,1}(:,j);nY(j)*Y_train{X_in==X_sub(i)}(:,j)]);max([Y_ns{i,1}(:,j);nY(j)*Y_train{X_in==X_sub(i)}(:,j)])])
            ind=ind+1;
            % saveas(gcf,['~/Desktop/y',num2str(j),'_n',num2str(i),'.png'])
        end
    end
end