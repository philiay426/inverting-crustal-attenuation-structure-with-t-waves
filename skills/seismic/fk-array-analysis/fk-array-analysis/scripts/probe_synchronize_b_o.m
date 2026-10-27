% probe_synchronize_b_o.m —— 切窗错位探针（Seizmo synchronize 会重定基准）
%
% 【为什么】`synchronize` 会把记录头重定基准（实测 b=-100,o=0 → b=0,o=100）。裸偏移 `cut(data,win)`
%   是相对【参考时刻】算的 ⇒ 凡输入 b ≠ 0 的脚本，实际切到的事件段比标称窗早 |b| 秒
%   （陆地计时输入 b = LAND_T0 = -100 ⇒ 早 100 s），症状 = 峰塌到近零慢度、方位反向，
%   而注入标定照样“通过”。跑本探针能一眼看出来。
%
% 【用法】改下面参数区 → conda run -n octave octave-cli --eval "run('<本文件绝对路径>')"
%   判读：`after sync` 的 o 非 0 且原文件 b ≠ 0 ⇒ 裸偏移脚本会错位 |b| 秒；
%         修法 = cut(data, o+win)（o 取 synchronize 之后的值）或裸偏移统一 +|b|。
% =============================================================================
cd('/home/sy/.hermes/workspace/seizmo/seizmo-0.6.16-Otgontenger'); startup_seizmo;

% ---- ★ 参数区（只改这里）----
INDIR   = '';                                   % ← 输入目录，例如 .../G1/04_landtime_fk/input_landtime
FILEPAT = '%s_land_%s.sac';                     % ← 文件名模板（STA, COMP）
STAS    = {'SA16','SA17','SA18','SA22','SA23','SA25','SA29'};
COMP    = 'T';
WIN     = [143.6 285.3];                        % ← 脚本里写的窗（裸偏移写法下的标称窗）
% -------------------------------

if isempty(INDIR), error('请先在参数区填 INDIR'); end
f = {}; for i = 1:numel(STAS), f{i} = fullfile(INDIR, sprintf(FILEPAT, STAS{i}, COMP)); end
d0 = readheader(f); d0 = readdata(d0);
b0 = getheader(d0(1),'b'); o0 = getheader(d0(1),'o'); dlt = getheader(d0(1),'delta');
fprintf('before sync : b=%.2f e=%.2f o=%.2f npts=%d\n', b0, getheader(d0(1),'e'), o0, numel(d0(1).dep));
d = synchronize(d0);
b1 = getheader(d(1),'b'); o1 = getheader(d(1),'o');
fprintf('after  sync : b=%.2f e=%.2f o=%.2f npts=%d\n', b1, getheader(d(1),'e'), o1, numel(d(1).dep));

% 逐台：裸偏移 WIN 切出的段（在重定后的坐标里）
for i = 1:numel(STAS)
    di = cut(d, i);
    v  = double(di.dep);
    b  = getheader(di,'b'); e = getheader(di,'e');
    fprintf('  %-6s cut(%.1f,%.1f) -> b=%.1f e=%.1f npts=%3d 段和=%.6g\n', ...
        STAS{i}, WIN(1), WIN(2), b, e, numel(v), sum(v));
    % 期望（若不加 o 修正）：应在参考相对 WIN(1)-b0 处切；这里报出“两者相差多少秒”
    if abs(b0) > 1e-6
        fprintf('         ⚠ 原文件 b=%.1f ≠ 0 ⇒ 裸偏移切到的段比“参考相对 = WIN−b_orig”早 %.1f s\n', b0, abs(b0));
        fprintf('           （相邻两台的段和不同只因起止样点，错位量对全台相同；用源记录扫等长窗可定位实际所在）\n');
        break;   % 错位结论与台站无关，报一台即可
    else
        fprintf('         ✓ 原文件 b=0：裸偏移语义与预期一致\n');
        break;
    end
end

% 修法示例（不改脚本时的手工验证）：把窗位整体 +|b0|
if abs(b0) > 1e-6
    dd = cut(d, WIN(1)+abs(b0), WIN(2)+abs(b0));
    fprintf('修正窗 +%.1f s 后：b=%.1f npts=%d 段和=%.6g（与上一行段和不同即“确实切在别处”）\n', ...
        abs(b0), getheader(dd,'b'), numel(dd.dep), sum(double(dd.dep)));
end
fprintf('判读：after sync 的 o=%.1f；裸偏移脚本一律改成 cut(data, o+win) 最稳。详见 references/land-time-fk.md §10\n', o1);
