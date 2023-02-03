% Skeleton code for simulation chain

% History:
%   2000-06-28  written /Stefan Parkvall
%   2001-10-22  modified /George Jongren
clear;
close all;

save_plots = true;
plot_ber = true;
plot_eye = true;
plot_cor = false;
window_size = 15; % Based on the correlation plots in sync, 15 seems best


% Initialization
EbN0_db = 0:10;                     % Eb/N0 values to simulate (in dB)
% EbN0_db = 10;                     
nr_bits_per_symbol = 2;             % Corresponds to k in the report
nr_guard_bits = 10;                 % Size of guard sequence (in nr bits)
                                    % Guard bits are appended to transmitted bits so
                                    % that the transients in the beginning and end
                                    % of received sequence do not affect the samples
                                    % which contain the training and data symbols.
nr_data_bits = 1000;                % Size of each data sequence (in nr bits)
nr_training_bits = 100;             % Size of training sequence (in nr bits)
nr_blocks = 50;                     % The number of blocks to simulate
Q = 8;                              % Number of samples per symbol in baseband

% Define the pulse-shape used in the transmitter. 
% Pick one of the pulse shapes below or experiemnt
% with a pulse of your own.
pulse_shape = ones(1, Q);
% pulse_shape = root_raised_cosine(Q);

% Matched filter impulse response. 
mf_pulse_shape = fliplr(pulse_shape);


% Loop over different values of Eb/No.
nr_errors = zeros(1, length(EbN0_db));   % Error counter
nr_errors_2_path = zeros(1, length(EbN0_db));

for snr_point = 1:length(EbN0_db)
  % Loop over several blocks to get sufficient statistics.
  for blk = 1:nr_blocks
    %%%
    %%% Transmitter
    %%%

    % Generate training sequence.
    b_train = training_sequence(nr_training_bits);
    
    % Generate random source data {0, 1}.
    b_data = random_data(nr_data_bits);

    % Generate guard sequence.
    b_guard = random_data(nr_guard_bits);
 
    % Multiplex training and data into one sequence.
    b = [b_guard b_train b_data b_guard];
    
    % Map bits into complex-valued QPSK symbols.
    d = qpsk(b);
  
    % Upsample the signal, apply pulse shaping.
    tx = upfirdn(d, pulse_shape, Q, 1);

    %%%
    %%% AWGN Channel
    %%%
    
    % Compute variance of complex noise according to report.
    sigma_sqr = norm(pulse_shape)^2 / nr_bits_per_symbol / 10^(EbN0_db(snr_point)/10);
    
    % Create noise vector.
    n = sqrt(sigma_sqr/2)*(randn(size(tx))+1j*randn(size(tx)));

    % Received signal.
    rx = tx + n;

    %%%
    %%% Receiver
    %%%
    
    % Matched filtering.
    mf=conv(mf_pulse_shape,rx);
    
    % Synchronization. The position and size of the search window
    % is here set arbitrarily. Note that you might need to change these
    % parameters. Use sensible values (hint: plot the correlation
    % function used for syncing)! 
    t_start=1+Q*nr_guard_bits/2;
    t_end=t_start + window_size;
    t_samp = sync(mf, b_train, Q, t_start, t_end, plot_cor && blk == 50);    
    % Down sampling. t_samp is the first sample, the remaining samples are all
    % separated by a factor of Q. Only training+data samples are kept.
    r = mf(t_samp:Q:t_samp+Q*(nr_training_bits+nr_data_bits)/2-1);

    % Phase estimation and correction.
    phihat = phase_estimation(r, b_train);
    r = r * exp(-1j*phihat);

    % Make decisions. Note that dhat will include training sequence bits
    % as well.
    bhat = detect(r);
    
    % Count errors. Note that only the data bits and not the training bits
    % are included in the comparison. The last data bits are missing as well
    % since the whole impulse response due to the last symbol is not
    % included in the simulation program above.
    temp=bhat(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
    nr_errors(snr_point) = nr_errors(snr_point) + sum(temp);



    %%%% ISI 2 Path

    tx_2_path = multipath(tx,8);
    
    % Compute variance of complex noise according to report.
    sigma_sqr_2_path = norm(pulse_shape)^2 / nr_bits_per_symbol / 10^(EbN0_db(snr_point)/10);
    
    % Create noise vector.
    n_2_path = sqrt(sigma_sqr_2_path/2)*(randn(size(tx_2_path))+1j*randn(size(tx_2_path)));

    % Received signal.
    rx_2_path = tx_2_path + n_2_path;

    %%%
    %%% Receiver
    %%%
    
    % Matched filtering.
    mf_2_path = conv(mf_pulse_shape,rx_2_path);
    
    % Synchronization. The position and size of the search window
    % is here set arbitrarily. Note that you might need to change these
    % parameters. Use sensible values (hint: plot the correlation
    % function used for syncing)! 
    t_start_2_path = 1+Q*nr_guard_bits/2;
    t_end_2_path = t_start_2_path + window_size;
    t_samp_2_path = sync(mf_2_path, b_train, Q, t_start_2_path, t_end_2_path, plot_cor && blk == 50);    
    % Down sampling. t_samp is the first sample, the remaining samples are all
    % separated by a factor of Q. Only training+data samples are kept.
    r_2_path = mf_2_path(t_samp_2_path:Q:t_samp_2_path+Q*(nr_training_bits+nr_data_bits)/2-1);

    % Phase estimation and correction.
    phihat_2_path = phase_estimation(r_2_path, b_train);
    r_2_path = r_2_path * exp(-1j*phihat_2_path);

    % Make decisions. Note that dhat will include training sequence bits
    % as well.
    bhat_2_path = detect(r_2_path);
    
    % Count errors. Note that only the data bits and not the training bits
    % are included in the comparison. The last data bits are missing as well
    % since the whole impulse response due to the last symbol is not
    % included in the simulation program above.
    temp_2_path = bhat_2_path(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
    nr_errors_2_path(snr_point) = nr_errors_2_path(snr_point) + sum(temp_2_path);

    % Next block.
  end
  % Next Eb/No value.
end

% Compute the BER. 
BER = nr_errors / nr_data_bits / nr_blocks;
BER_2_path = nr_errors_2_path / nr_data_bits / nr_blocks;

if plot_ber
    figure
    semilogy(EbN0_db, BER, "b")
    hold on
    semilogy(EbN0_db, BER_2_path, "r")
    legend("AWGN", "2 Path")
    xlabel('Eb/N0');
    ylabel('Pe');
    if save_plots
        saveas(gcf, "./images/BER_2_path.png");
    end

end

if plot_eye
    eyediagram(r,2);
    if save_plots
        saveas(gcf, "./images/eye.png");
    end
    eyediagram(r_2_path,2);
    if save_plots
        saveas(gcf, "./images/eye_2_path.png");
    end
   
end
