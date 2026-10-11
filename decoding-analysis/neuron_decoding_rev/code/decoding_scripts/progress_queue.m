function q = progress_queue(label, total)
%{
Filename: progress_queue.m
Author:   Zitong Wang
Date:     2026-10-11

Description:
    Text progress for a parfor loop. Returns a parallel.pool.DataQueue; call
    send(q, 1) once per finished iteration inside the parfor. Roughly every
    10% of `total` (and at the end) one line is printed on the client:
        <label>: 36/120 calls (30%) | 4.1 min elapsed, ~9.5 min left
    Plain text rather than a waitbar, so it also works in batch runs and is
    kept in the master script's log file.
%}

q = parallel.pool.DataQueue;
n = 0;
t0 = tic;
step = max(1, ceil(total / 10));
next = step;
afterEach(q, @tick);

    function tick(~)
        n = n + 1;
        if n >= next || n == total
            el = toc(t0);
            fprintf('    %s: %d/%d calls (%3.0f%%) | %.1f min elapsed, ~%.1f min left\n', ...
                label, n, total, 100 * n / total, el / 60, el / n * (total - n) / 60);
            next = next + step;
        end
    end
end
