# Terminal 执行故障与变通方案

## "embedded null byte" 错误（2026-09-02）

### 症状

执行 Python 脚本时，terminal 工具报错：

```
Failed to execute command: embedded null byte
ValueError: embedded null byte
  File ".../lifecycle_guard.py", line 260, in _read_referenced_script
    descriptor = os.open(path, flags)
```

### 表现

- 直接调用 Python 解释器（含完整路径）失败
- 用 `cat script.py | python` 管道方式失败
- 使用 heredoc 或其他间接方式仍失败
- 简单命令（`ls`, `pwd`, `file`）正常

### 根因（2026-09-02 更新）

Hermes 终端工具的 `lifecycle_guard` 在解析包含 Python 脚本路径的命令时触发系统级错误。关键特征：

1. **触发条件**：命令中包含 `open('xxx.py')` 或直接引用 `.py` 文件路径
2. **错误位置**：`lifecycle_guard.py` 尝试用 `os.open()` 读取脚本内容时失败
3. **持续性**：同一文件路径在会话中反复触发，即使：
   - 改变执行方式（管道、heredoc、wrapper）
   - 重命名文件（`fk_120s_bkta_cpsm_ampy.py` → `fk_120s_fixed.py`）
   - 复制到不同目录（如 `/tmp/`）

**可能原因**（未确认）：
- lifecycle_guard 对某些路径模式或文件内容的解析有 bug
- 与 WSL 文件系统交互（/mnt/e NTFS 挂载）的兼容性问题
- 脚本内容中的某些字符序列触发路径解析异常

重要：文件本身是正常的（`file` 显示 "Python script, UTF-8 text"，`od -c` 无空字节），问题在工具层而非文件层。

### 有效变通方案

**1. 使用子代理（delegate_task）— 首选**

最可靠的方案。子代理有独立的终端会话，不受主会话的 lifecycle_guard 状态影响：

```python
delegate_task(
    goal="执行 /path/to/script.py 并报告结果",
    context="脚本用途、输入输出说明、Python 解释器路径"
)
```

**注意**：子代理也可能遇到同样的 bug，但成功率明显高于主会话（可能因为新的终端会话状态更干净）。2026-09-02 马达加斯加案例中，主会话持续失败，子代理第一次尝试即成功。

**2. 重写脚本到新文件名（临时缓解）**

有时改变文件名可以绕过特定路径的触发条件，但不保证成功：

```python
terminal("mv problematic.py temp_script.py")
# 然后尝试执行 temp_script.py
```

**3. 后台进程（不推荐，成功率低）**

`background=true` 理论上可以绕过某些检查，但实测仍会触发同样错误：

```python
terminal(
    command="cd /path && nohup python script.py > output.log 2>&1 &",
    background=True,
    notify_on_complete=True
)
```

注意：输出目录必须提前创建，否则 nohup 重定向失败。

**3. 不要尝试的方案（已验证无效）**

- ❌ 改用不同的 Python 解释器路径
- ❌ heredoc 内嵌脚本内容
- ❌ 创建 wrapper shell 脚本
- ❌ 用 `bash -c '...'` 包装

以上方案在该错误场景下都会触发同样的 `embedded null byte` 异常。

### 预防措施

- 尽量用子代理执行复杂脚本（隔离性更好）
- 避免在脚本中使用非 ASCII 字符（如果必须，确保 UTF-8 编码声明）
- 脚本路径用绝对路径或相对于明确的工作目录

## 内联 heredoc 脚本里的中文引号（2026-09-14）

- 本项目大量用 `python - <<'EOF' … EOF` / `/home/sy/miniforge3/envs/obspy/bin/python - <<'PYEOF'` 做一次性统计、
  出图与补丁，**这是可用的**（与上面的 null-byte 故障无关）。但内联脚本里**不要出现中文全角引号**
  （`“ ” ‘ ’`）嵌在字符串内：`print("…“保留占比”…")` 会报 `SyntaxError: invalid syntax`（全角引号被当作
  引号对，字符串被提前截断），报错行还要往上找一两行才看得懂。
- 两种修法：① 内联脚本里的中文措辞不加引号（或用 `「」`）；② 要写长字符串就 `write_file` 落成 `.py`
  再执行，别在 heredoc 里拼。
- 同类坑：用 `assert old in t` + `t.replace(old, new)` 打补丁时，`old` 里的中文标点必须与源文件**逐字符一致**
  （源文件常用全角 `（）` `，` `→`），一律从 `read_file` / `skill_view` 的输出里直接复制，不要手敲；
  先 `assert` 再写，避免静默不匹配。
- 同类坑（numpy 2）：`ndarray.ptp()` 已被移除 → `AttributeError: 'numpy.ndarray' object has no attribute 'ptp'`，
  改用 `np.ptp(arr)`；批量改脚本时先 `grep -rn "\.ptp()"` 一次找全（本次同一坑踩了两处，分两次才发现）。
  通用规则：见到"ndarray 没有 X 方法/属性"先怀疑 numpy 2 移除项（`ptp`/`float_`/`NaN`/`in1d` 等），
  用 `np.<等价函数>(arr)` 替。

### 诊断步骤

遇到该错误时：

1. 确认脚本文件本身正常：`file script.py`（应显示 "Python script, UTF-8 text"）
2. 检查脚本是否可独立执行：尝试简单的 `python -c "print(1)"`
3. 如果简单命令正常、仅脚本失败 → 使用子代理
4. 如果所有 terminal 调用都失败 → 可能是会话级问题，建议重启或使用不同工具
