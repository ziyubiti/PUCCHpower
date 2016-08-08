
function hHARQACKResultsV2(SNRdB, pFalse, pMissed)

    figure;
    % Probability of false detection
    semilogy(SNRdB, pFalse, 'k-*', 'LineWidth', 2, 'MarkerSize', 7); 
    grid on;
    hold on;
    % Probability of missed detection
    semilogy(SNRdB, pMissed, 'b-o', 'LineWidth', 2, 'MarkerSize', 7); 
    hold on;
    % Plot target probability
%     plot(13.8, 0.01, 'rx', 'Linewidth', 2, 'MarkerSize', 7);
    xlabel('SNR (dB)');
    ylabel('Probability of false/missed detection');
    title('ACK/NACK on PUSCH ');   
    axis([SNRdB(1)-0.1 SNRdB(end)+0.1 0.001 1.0]);
    legend('False probability', ...
        'Miss probability', 'target');

end