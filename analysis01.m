%%
fs = 16000;

%%{
f_lb = 60;
f_ub = 250;
dataname = 'Yuto_Ozaki_Japanese_Japanese_Traditional_Asatoya-Yunta_20220209_desc';
%}

%{
f_lb = 70;
f_ub = 400;
dataname = 'Yuto_Ozaki_Japanese_Japanese_Traditional_Asatoya-Yunta_20220209_song';
%}

audiofile = strcat(dataname, '.wav');
mkdir(strcat('./fig/', dataname, '/'));

%%
[s, fs_orig] = audioread(audiofile);
if size(s, 2) == 2
    s = mean(s, 2);
end

%%
s = h_downsampling(s, fs_orig, fs);

%%
addpath(strcat(userpath, filesep, 'lib2', filesep, 'ewtbpfilt', filesep));
gpumode = true;
y = gather(ewtbandpass(s, f_lb, f_ub, fs, gpumode));

%%
segment = h_voicingknn(y, round(fs*0.01), round(fs*0.02));

%%
be = 30;
gam = 9;
frange = [25, 3200];
voice = 120;
J = 3;
N = 20;
L_integer = round(h_lengthcheck(be, gam, frange(1), fs)*1.1);
offset = ceil(L_integer/2);

oct = log2(frange(end)/frange(1));
F = frange(1).*2.^((0:(oct*voice))'./voice);
q = 8;

for k=75:size(segment, 1)
    % ConceFT
    z = s(segment(k, 1) - offset:segment(k, 2) + offset);
    [W, T, E_Omg] = h_conceFT(z', be, gam, frange, voice, fs, J, N);
    W = W(:, 1 + offset:end - offset);
    T = T(:, 1 + offset:end - offset);
    E_Omg = E_Omg(:, 1 + offset:end - offset);
    
    % De-shape
    U = h_gausscorr(abs(T), voice);
    U = U./max(abs(U(:)));
    D = abs(T).*2.*exp(q.*U);
    
    % Viterbi
    V = D;
    V(D == 0) = -Inf;
    R = h_viterbi(V, 1, Inf);
    IF = R.*0;
    for n=1:numel(IF)
        IF(n) = E_Omg(R(n), n);
    end

    % plot
    fobj = figure(1);
    clf; cla;
    fobj.Position = [40, 580, 1600, 400];

    subplot(1, 4, 1);
    p = pcolor(1:size(T, 2), F, log(abs(T).^2)); set(gca, 'YDir', 'normal');
    p.EdgeColor = 'none';
    title(['k = ', num2str(k, '%d'), ', \beta = ', num2str(be, '%3.3f'), ', \gamma = ', num2str(gam, '%3.3f'),...
        ', J = ', num2str(J, '%d'), ', N = ', num2str(N, '%d')]);
    ylim([25, 1000]);

    subplot(1, 4, 2);
    p = pcolor(1:size(T, 2), F, log(abs(W).^2)); set(gca, 'YDir', 'normal');
    p.EdgeColor = 'none';
    ylim([25, 1000]);

    subplot(1, 4, 3);
    p = pcolor(1:size(T, 2), F, U); set(gca, 'YDir', 'normal');
    p.EdgeColor = 'none';
    ylim([25, 1000]);

    subplot(1, 4, 4);
    p = pcolor(1:size(T, 2), F, log(D)); set(gca, 'YDir', 'normal');
    p.EdgeColor = 'none';
    ylim([25, 1000]);
    hold on
    plot(1:size(T, 2), F(R), 'Color', 'R');
    hold off
    
    drawnow

    saveas(fobj, strcat('./fig/', dataname, '/', dataname, ' - ', num2str(k, '%d'), '.jpg'));
end