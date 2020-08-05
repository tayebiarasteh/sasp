% load signal
disp('Reading input signal...')
[signal, fs] = wavread('noisySignal');
signal = resample(signal, 16e3, fs);

% apply center clipper
disp('Applying center clipper (threshold eta=0.01)...')
filteredSignal = cclipper(signal, 0.01);

% listen to the input signal
disp('Listening to the input signal...')
soundsc(signal, 16e3)

% listen to the filtered signal
disp('Listening to the output signal...')
soundsc(filteredSignal, 16e3)
