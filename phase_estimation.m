function phihat = phase_estimation(r, b_train)
% phihat = phase_estimation(r, b_train)
%
% Phase estimator using the training sequence. The phase estimate is
% obtained by minimizing the norm of the difference between the known
% transmitted QPSK-modulated training sequence and the received training
% part. NB! There are other ways of estimating the phase, this is just
% one example.
%
% Input:
%   r       = received baseband signal
%   b_train = the training sequence bits
%
% Output:
%   phihat     = estimated phase
phihat = 0;
qpsk_train = qpsk(b_train);
n = r(1:length(qpsk_train));
min_g = norm(n-qpsk_train);
% Define the range and the step of sequence
for phase = -pi:0.1:pi
    r_sym = n*exp(-1j*phase);
    min_l = norm(r_sym-qpsk_train);
    if min_l < min_g
        min_g = min_l;
        phihat = phase;
    end
end

