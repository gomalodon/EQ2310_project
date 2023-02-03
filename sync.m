function t_samp = sync(mf, b_train, Q, t_start, t_end)
% t_samp = sync(mf, b_train, Q, t_start, t_end)
%
% Determines when to sample the matched filter outputs. The synchronization
% algorithm is based on cross-correlating a replica of the (known)
% transmitted training sequence with the output from the matched filter
% (before sampling). Different starting points in the matched filter output
% are tried and the shift yielding the largest value of the absolute value
% of the cross-correlation is chosen as the sampling instant for the first
% symbol.
%
% Input:
%   mf            = matched filter output, complex baseband, I+jQ
%   b_train       = the training sequence bits
%   Q             = number of samples per symbol
%   t_start       = start of search window
%   t_end         = end of search window
%
% Output:
%   t_samp = sampling instant for first symbol
L = t_end - t_start;
c = qpsk(b_train); 
correlations = zeros(1,L);
for i=t_start:t_end
    temp = mf(i+((0:length(c)-1).*Q));
    correlations(i - t_start + 1) = abs(temp*c');
end
[top,index] = max(correlations);
t_samp = index + t_start - 1;

if false
    stem(t_start:t_end,correlations, "b")
    hold on
    stem(t_samp,top, "r")
end
