clear all; close all; clc

Responses = [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1];
MaxResponse = 1;
BackgroundProb = .5;

runanalysisv2(Responses, MaxResponse, BackgroundProb)

load resultsindividual.mat

fontsize1 = 20;
linewidth1 = 1;

 t=1:size(p,2)-1; 

 figure(1);  clf;

 %plot learning curve
 subplot(111);  
 h = plot(t, pmid(2:end),'b-'); set(h, 'LineWidth',2);
 hold on;
 plot(t, p05(2:end),'b', t, p95(2:end), 'b', 'LineWidth', linewidth1);
 if(MaxResponse == 1)
      hold on; [y, x] = find(Responses > 0);
      h = plot(x,y,'o'); set(h, 'MarkerFaceColor',[.9 .9 .9]);
      set(h, 'MarkerEdgeColor', 'k' ,'MarkerSize', 4);
      hold on; [y, x] = find(Responses == 0);
      h = plot(x,zeros(size(x)),'o'); set(h, 'MarkerFaceColor', [.9 .9 .9]);
      set(h, 'MarkerEdgeColor', 'k','MarkerSize', 4);
      axis([1 t(end)  0 1.0]);
 else
      hold on; plot(t, Responses./MaxResponse,'ko');
      axis([1 t(end)  0 1]);
 end
 plot([1 t(end)], [BackgroundProb  BackgroundProb ], 'b', 'LineWidth', linewidth1);
 title(['Learning trial = ' num2str(cback)],'fontsize',fontsize1);
 xlabel('Trial Number','fontsize',fontsize1)
 ylabel('Probability of a correct response','fontsize',fontsize1)
 set(gca,'tickdir','out','fontsize',fontsize1), box off

