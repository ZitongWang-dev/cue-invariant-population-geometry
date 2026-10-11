function m = transfer_metrics(self, pt, ctrl, nopt)
%{
Filename: transfer_metrics.m
Author:   Zitong Wang
Date:     2026-10-11

Description:
    Pooled summary of transfer-decoding accuracies, as reported in the rebuttal
    and the simulation write-ups: every quantity is a ratio of POOLED MEANS
    (mean PT / mean control), never a mean of per-repeat ratios. Inputs can have
    any shape (e.g. pairs x draws x folds); NaN entries are ignored.

Inputs:
    self, pt, ctrl : self-decoding, PT (rotation only) and corrected-control accuracies
    nopt           : optional non-PT accuracies

Output: struct m with fields
    self, pt, ctrl, nopt   pooled means
    ratio                  pt / ctrl
    diff                   pt - ctrl
    headroom               (pt - ctrl) / (self - ctrl): the share of the
                           above-floor decodability that transfers
%}

m.self = mean(self(:), 'omitnan');
m.pt   = mean(pt(:),   'omitnan');
m.ctrl = mean(ctrl(:), 'omitnan');
if nargin > 3 && ~isempty(nopt)
    m.nopt = mean(nopt(:), 'omitnan');
else
    m.nopt = NaN;
end
m.ratio    = m.pt / m.ctrl;
m.diff     = m.pt - m.ctrl;
m.headroom = (m.pt - m.ctrl) / (m.self - m.ctrl);
end
