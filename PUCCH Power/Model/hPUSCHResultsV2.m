

function hPUSCHResultsV2(SNRIn, NFrames, trBlkSizes, throughput, bitThroughput)

    figure;
    plot(SNRIn, mean(bitThroughput, 2)/1000,'-b*');
    title(['Throughput for ', num2str(NFrames) ' Frame(s)']);
    xlabel('SNR (dB)'); ylabel('Throughput (Mbps)');
    grid on;
    hold on;
    plot(SNRIn, mean(trBlkSizes)*0.95*ones(1, numel(SNRIn))/1000, '--rs',SNRIn, mean(trBlkSizes)*0.7*ones(1, numel(SNRIn))/1000, '--cs',SNRIn, mean(trBlkSizes)*0.3*ones(1, numel(SNRIn))/1000, '--ms');
    set(gca, 'YLim', [0 mean(trBlkSizes)*1.05/1000]);
    legend('Simulation Result', '95 Percent Throughput','70 Percent Throughput', '30 Percent Throughput', ...
        'Location', 'SouthEast')
    
    figure;
    plot(SNRIn, throughput,'-b*');
    title(['Throughput for ', num2str(NFrames) ' Frame(s)']);
    xlabel('SNR (dB)'); ylabel('Throughput (%)');
    grid on;
    hold on;
    plot(SNRIn, 95*ones(1, numel(SNRIn)), '--rs',SNRIn, 70*ones(1, numel(SNRIn)), '--cs',SNRIn, 30*ones(1, numel(SNRIn)), '--ms');
    set(gca, 'YLim', [0 105]);
    legend('Simulation Result', '95 Percent Throughput','70 Percent Throughput','30 Percent Throughput', ...
        'Location', 'SouthEast');
    
    
end