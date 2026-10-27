% =============================================================================
% step3_fk_driver.m —— 步3 F-K 最小驱动（Octave + SEIZMO fkmap）
% -----------------------------------------------------------------------------
% 对应 skill：seismic/twave-crustal-attenuation §5 步 3；引擎实现见 skill
%   fk-array-analysis（fk4d/fkmap、参数公式、子阵选择判据都在那里）。
% 处理链：N sps BHZ 2-8 Hz 波形 -> synchronize -> interpolate 窗内
%   -> normalize -> |Hilbert| 包络 -> getspline(0.1) -> removemean -> taper
%   -> fkmap(smax, spts, frng)
% 台站 = FK_DATA 目录里全部 SAC（本脚本不认台名；子阵由目录内容决定）。
% 台太少（<4-5）会在上界频率以上混叠、峰跑到错方位 —— 判据见 step3_fk_demo.py 的打印。
%
% 运行（必须经 conda，直接调 octave-cli 会报 fileparts undefined）：
%   FK_SEIZMO=<seizmo 目录> FK_DATA=<本步>/data/g1 FK_TAG=sub1 T_OUTDIR=<本步> \
%   conda run -n octave octave-cli step3_fk_driver.m
% =============================================================================
prog = @(s) fprintf(stderr, '%s\n', s);
pkg load signal; pkg load splines;
SEIZMO_DIR = getenv('FK_SEIZMO');
if isempty(SEIZMO_DIR)
  error('请设置环境变量 FK_SEIZMO 指向 SEIZMO 安装目录');
end
cd(SEIZMO_DIR); startup_seizmo;

% ---- ★ 参数区（全部可用环境变量覆盖）----
TAG  = getenv('FK_TAG');  if isempty(TAG), TAG = 'sub1'; end
DATA = getenv('FK_DATA');
if isempty(DATA), error('请设置 FK_DATA 指向该子阵的 SAC 目录'); end
OUT  = getenv('T_OUTDIR'); if isempty(OUT), OUT = pwd; end
% 窗：必填。不同事件/台阵完全不同，不许有默认值
_w = str2num(getenv('FK_TWIN'));
if isempty(_w) || numel(_w) ~= 2
  error('请设置 FK_TWIN="W0,W1"（同步坐标下的分析窗，单位 s；窗长 = W1-W0）');
end
TWIN = _w; WINLEN = TWIN(2) - TWIN(1);
% 网格：优先 FK_SMAX；否则由预测慢度现算 smax = 1.4 * u_pred, u_pred = 111.195/FK_VT
SMAX = str2double(getenv('FK_SMAX'));
if isnan(SMAX)
  VT = str2double(getenv('FK_VT'));
  if isnan(VT)
    error('请设置 FK_VT（T 波视速度 km/s）或 FK_SMAX（网格半径 s/deg）；推荐先跑 fk_geometry.py');
  end
  SMAX = 1.4 * 111.195 / VT;
end
SPTS = str2double(getenv('FK_SPTS')); if isnan(SPTS), SPTS = 301; end      % 奇数，0 在格点
% 频带：默认 [1/窗长, min(FK_FMAX, C/smax)]；C = 111.195/(2*d_min) 由 fk_geometry.py 给出
_f   = str2num(getenv('FK_FRNG'));
FMAX = str2double(getenv('FK_FMAX')); if isnan(FMAX), FMAX = 0.04; end
CFK  = str2double(getenv('FK_C'));
if isempty(_f)
  if isnan(CFK), FRNG = [1/WINLEN FMAX]; else FRNG = [1/WINLEN min(FMAX, CFK/SMAX)]; end
else
  FRNG = _f;
end
prog(sprintf('窗 [%g %g] s 长 %g s | smax=%g spts=%d frng=[%.5f %.5f] Hz', ...
     TWIN(1), TWIN(2), WINLEN, SMAX, SPTS, FRNG(1), FRNG(2)));
SPS  = str2double(getenv('FK_SPS')); if isnan(SPS), SPS = 20; end         % 重采样率
if ~exist(OUT, 'dir'); mkdir(OUT); end

% ---- 1) 读入该子阵 2-8 Hz 波形 ----
sf = readseizmo(fullfile(DATA, '*.SAC*'));
sf = sortbyfield(sf, 'kstnm');
stn = getheader(sf, 'kstnm');
prog(['子阵: ' strjoin(cellstr(stn), ' ') '  (' num2str(numel(sf)) ' 台)']);

% ---- 2) 时间对齐 + 切窗 ----
ssfb = synchronize(sf);
ssfb = interpolate(ssfb, SPS, 'spline', TWIN(1), TWIN(2));
[npts, delta] = getheader(ssfb(1), 'npts', 'delta');
prog(sprintf('窗内 npts=%d delta=%g (nspts=%d, df=%g Hz)', ...
     npts, delta, 2^nextpow2(npts), 1/(2^nextpow2(npts)*delta)));

% ---- 3) 预处理后三步 + 包络（必须取包络）----
ssfb = normalize(ssfb);
ssfb = envelope(ssfb);
ssfb = getspline(ssfb, 0.1);
ssfb = removemean(ssfb);
ssfb = taper(ssfb);

% ---- 4) 波束成形 ----
samp = fkmap(ssfb, SMAX, SPTS, FRNG);
x = samp.x(:)'; y = samp.y(:); B = samp.beam;

[~, idx] = max(B(:)); [iy, ix] = ind2sub(size(B), idx);
prog(sprintf('峰 |u| = %.2f s/deg   BAZ = %.1f deg   (格点 x=%d y=%d)', ...
     hypot(x(ix), y(iy)), mod(atan2d(x(ix), y(iy)), 360), ix, iy));

% ---- 5) 存网格（Python 侧出图）----
fid = fopen(fullfile(OUT, ['beam_' TAG '.txt']), 'w');
fprintf(fid, '# %s %d 台子阵 | 2-8 Hz BHZ 包络 + getspline(0.1) | smax=%g spts=%d frng=[%g %g]\n', TAG, numel(sf), SMAX, SPTS, FRNG(1), FRNG(2));
fprintf(fid, '# 列: x=E/W 慢度(s/deg)  y=N/S 慢度(s/deg)  beam(dB, 峰归一到 0)\n');
for i = 1:numel(y)
  for j = 1:numel(x)
    fprintf(fid, '%12.6f %12.6f %12.6f\n', x(j), y(i), B(i, j));
  end
end
fclose(fid);

% ---- 6) 存送进 fkmap 的波形（画"时窗落在波形上"用）----
fid2 = fopen(fullfile(OUT, ['fkinput_' TAG '.txt']), 'w');
fprintf(fid2, '# %s fkmap 输入: 2-8 Hz 包络 -> getspline(0.1) -> removemean -> taper\n', TAG);
fprintf(fid2, '# nsta npts delta t_start\n');
fprintf(fid2, '# %d %d %g %g\n', numel(ssfb), npts, delta, TWIN(1));
for i = 1:numel(ssfb)
  nm = getheader(ssfb(i), 'kstnm'); nm = nm{1};
  fprintf(fid2, '%s', nm);
  fprintf(fid2, ' %.6f', ssfb(i).dep);
  fprintf(fid2, '\n');
end
fclose(fid2);
prog('FK_DONE');
