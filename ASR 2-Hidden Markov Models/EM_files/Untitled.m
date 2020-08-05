h1 = waitbar(0, 'Waitbar 1', 'Units', 'normalized', 'Position', [0.25 0.4 0 0]);
for k = 1:10
    waitbar(k/10, h1);
    pause(0.1)
    h2 = waitbar(0, 'Waitbar 2', 'Units', 'normalized', 'Position', [0.5 0.4 0 0]);
    for m = 1:20
        waitbar(m/20, h2);
        pause(0.1)
    end
    close(h2)
end
