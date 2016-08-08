

function hPRACHDetectionResultsV2(SNRdB, numSubframes, P)

    figure;
    plot(SNRdB,P,'b-o','LineWidth',2,'MarkerSize',7);
    title(['Detection Probability for ', ... 
        num2str(numSubframes) ' subframe(s)'] );
    xlabel('SNRdB'); ylabel('Detection Probability');
    grid on;
    hold on;
%     plot(-8.0,0.99,'rx','LineWidth',2,'MarkerSize',7);
    legend('PRACH Miss Performence');
%         'Target 99% Probability',...
%         'Location','SouthEast');
    axis([SNRdB(1)-0.1 SNRdB(end)+0.1 0.0 1.05]) 

end