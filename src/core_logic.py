# ==============================================================================
# 主程序 (EIS Analysis & Fitting)
# ==============================================================================
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import os
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

from impedance.models.circuits import CustomCircuit
from scipy.signal import savgol_filter
from scipy.stats import linregress, f as f_dist 
from matplotlib.gridspec import GridSpec

plt.rcParams['font.sans-serif'] = ['Arial', 'SimHei', 'Microsoft YaHei']


plt.rcParams['pdf.fonttype'] = 42  # AI可编辑字体
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.size'] = 12     # 统一全局字号
# --- 主程序配置区 ---
plt.rcParams['font.sans-serif'] = ['Arial', 'SimHei', 'Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False



# --- 主程序函数定义 ---

def calculate_f_test_pvalue(z_obs, z_pred_simple, z_pred_complex, p_simple=3, p_complex=6):
    """
    执行嵌套模型的纯粹 F-检验 (引入模值权重 Modulus Weighting)。
    使用相对百分比误差，防止低频尾巴掩盖高频微小半圆。
    """
    n = len(z_obs)
    if n <= p_complex: return 1.0 
    
    # 计算相对残差平方和
    rss_simple = np.sum((np.abs(z_obs - z_pred_simple) / np.abs(z_obs))**2)
    rss_complex = np.sum((np.abs(z_obs - z_pred_complex) / np.abs(z_obs))**2)
    
    if rss_simple <= rss_complex: return 1.0 
    
    df1 = p_complex - p_simple  
    df2 = n - p_complex         
    
    F_stat = ((rss_simple - rss_complex) / df1) / (rss_complex / df2)
    return f_dist.sf(F_stat, df1, df2)

def manual_kk_test(freq, z):
    """ 动态 KK 校验 (带物理边界约束) """
    try:
        # 安全处理频率范围
        f_min, f_max = np.max([np.min(freq), 1e-6]), np.max(freq)
        decades = max(1, int(np.log10(f_max) - np.log10(f_min)))
        n_rc = min(decades, 6) # 自适应 RC 数量
        
        circuit_str = 'R0' + ''.join([f'-p(R{i},C{i})' for i in range(1, n_rc + 1)])
        f_centers = np.logspace(np.log10(f_min), np.log10(f_max), n_rc)
        
        r_guess = np.abs(z[0]) / n_rc  
        initial_guess = [max(np.real(z[0]), 1e-6)] 
        for fc in f_centers:
            tau = 1.0 / (2 * np.pi * fc)
            initial_guess.extend([r_guess, tau / r_guess])
            
        kk_circuit = CustomCircuit(circuit_str, initial_guess=initial_guess)
        weights = np.concatenate((np.abs(z), np.abs(z)))
        
        # 施加严格的非负物理边界，避免矩阵奇异崩溃
        bounds_low = [1e-15] * len(initial_guess)
        bounds_high = [np.inf] * len(initial_guess)
        
        kk_circuit.fit(freq, z, sigma=weights, bounds=(bounds_low, bounds_high))
        
        z_kk_full = kk_circuit.predict(freq)
        res_real = (np.real(z - z_kk_full)) / np.abs(z)
        res_imag = (np.imag(z - z_kk_full)) / np.abs(z)
        return res_real, res_imag, True
    except:
        return np.ones_like(freq), np.ones_like(freq), False
def load_and_clean_data(filepath):
    try:
        header_row = 0
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for i, line in enumerate(f):
                if any(k in line for k in ["Freq", "Frequency"]) and any(k in line for k in ["Z'", "Re(Z)"]):
                    header_row = i
                    break
        df = pd.read_csv(filepath, skiprows=header_row, sep=',')
        df.columns = [c.strip() for c in df.columns]
        col_f = [c for c in df.columns if "Freq" in c or "Hz" in c][0]
        col_zr = [c for c in df.columns if "Z'" in c or "Re" in c][0]
        col_zi = [c for c in df.columns if 'Z"' in c or "Im" in c][0]
        f_raw, z_raw = df[col_f].values, df[col_zr].values + 1j * df[col_zi].values
        sort_idx = np.argsort(f_raw)[::-1]
        freq, z = f_raw[sort_idx], z_raw[sort_idx]
        mask = (np.real(z) > 0) & (np.imag(z) <= 0) & np.isfinite(np.real(z)) & np.isfinite(np.imag(z))
        return freq[mask], z[mask], freq[~mask], z[~mask]
    except: 
        return None, None, None, None

def filter_bad_z_prime(z):
    if z is None or len(z) == 0: return None, []
    zr = np.real(z)
    n = len(zr)
    valid_indices = []
    bad_indices = []
    search_limit = min(20, int(n * 0.2) + 1)
    start_idx = np.argmin(zr[:search_limit]) if search_limit > 0 else 0
    for k in range(start_idx): bad_indices.append(k)
    current_max_zr = -1e9
    if start_idx < n:
        current_max_zr = zr[start_idx]
        valid_indices.append(start_idx)
    for i in range(start_idx + 1, n):
        val = zr[i]
 # [安全加固]：取高频前10个点（如果不够就有几个取几个）的实部极差
        hf_diffs = np.abs(np.diff(zr[:min(10, n)]))
        # 用中位数估计本底抖动，乘以3作为基础容差，保底 1e-6 防止全0
        base_noise = np.median(hf_diffs) * 3.0 if len(hf_diffs) > 0 else 1e-6
        
        # 将原先的 2.0 替换为 base_noise
        tolerance = max(base_noise, current_max_zr * 0.03)
        if val >= current_max_zr - tolerance:
            valid_indices.append(i)
            if val > current_max_zr: current_max_zr = val
        else:
            bad_indices.append(i)
    return z[valid_indices], bad_indices

def check_segment_quality(segment, min_len=5, min_amp=0.5):
    if len(segment) < min_len: return False, "太短"
    total_change = np.abs(segment[-1] - segment[0])
    if total_change < min_amp: return False, "幅度太小"
    diffs = np.diff(segment)
    is_rising = segment[-1] > segment[0]
    monotonic_steps = np.sum(diffs > 0) if is_rising else np.sum(diffs < 0)
    if monotonic_steps / len(diffs) < 0.75: return False, "非单调"
    if np.max(np.abs(diffs)) > total_change * 0.8: return False, "突变"
    return True, "合格"

def detect_shape(z_raw):
    z_clean, _ = filter_bad_z_prime(z_raw)
    if z_clean is None or len(z_clean) < 10: return "tail_only", 0
    y = -np.imag(z_clean)
    n = len(y)
    window = max(11, int(len(y)/8))
    if window % 2 == 0: window += 1
    if window >= n: window = max(3, n // 2 - 1) if n > 3 else 3
    try: y_smooth = savgol_filter(y, window, 2)
    except: y_smooth = y
    
    total_range = np.max(y_smooth) - np.min(y_smooth)
    search_limit = int(n * 0.9)
    candidates = []
    
    for i in range(2, search_limit):
        if y_smooth[i] <= y_smooth[i-1] and y_smooth[i] <= y_smooth[i+1]:
            valley_idx = i
            peak_idx = np.argmax(y_smooth[:valley_idx])
            segment_rise = y_smooth[0 : peak_idx + 1]
            segment_fall = y_smooth[peak_idx : valley_idx + 1]
            rise_ok, _ = check_segment_quality(segment_rise)
            fall_ok, _ = check_segment_quality(segment_fall)
            
            if rise_ok or fall_ok:
                tail_rise = np.max(y_smooth[valley_idx:]) - y_smooth[valley_idx]
                if tail_rise > max(total_range * 0.05, 1.0):
                    
                    n_rise = len(segment_rise) if rise_ok else 0
                    n_fall = len(segment_fall) if fall_ok else 0
                    
                    # ========================================================
                    # 自适应特征权重函数
                    # 采用全局有效采样点数 (n) 
                    # 与局部有效分布 (n_rise + n_fall) 动态计算权重。
                    # ========================================================
                    # 动态权重函数 f(N_total, N_local)：
                    # 使得归一化的拓扑高度 (0~1之间) 能够动态映射到与实际采样点数同等量级的评价尺度上。
                    lambda_weight = (n_rise + n_fall) * 0.5 + n * 0.2
                    
                    # 除以 total_range 实现无量纲化 (加 1e-9 防止除零报错)
                    normalized_z_drop = (y_smooth[peak_idx] - y_smooth[valley_idx]) / (total_range + 1e-9)
                    
                    # score 具备对不同设备采样密度与阻抗量级的自适应泛化能力
                    score = n_rise + n_fall + (lambda_weight * normalized_z_drop)
                    
                    candidates.append({'v_idx': valley_idx, 'score': score})
                    
    if candidates:
        best = sorted(candidates, key=lambda x: x['score'], reverse=True)[0]
        return "semicircle_with_tail", best['v_idx']
    else:
        return "tail_only", 0
def geometric_tail_analysis(z_tail):
    x, y = np.real(z_tail), -np.imag(z_tail)
    if len(x) < 4: return np.min(x), 0.7

    window_size = 4
    best_slope = None
    best_intercept = None
    found = False

    search_limit = max(window_size, int(len(x) * 0.6))
    for i in range(search_limit - window_size + 1):
        x_sub, y_sub = x[i : i + window_size], y[i : i + window_size]
        slope, intercept, r_val, _, _ = linregress(x_sub, y_sub)
        if r_val**2 > 0.95 and slope > 0.1:
            best_slope, best_intercept = slope, intercept
            found = True
            break  

    if not found:
        slope, intercept, _, _, _ = linregress(x[:int(len(x)*0.5)+1], y[:int(len(x)*0.5)+1])
        best_slope = slope if slope > 0 else 1.0
        best_intercept = intercept

    angle = np.arctan(best_slope)
    n_est = angle / (np.pi / 2)
    n_est = max(0.15, min(1.05, n_est))

    if best_slope > 2:
        x_intercept = np.min(x)
    else:
        math_intercept = -best_intercept / best_slope
        x_intercept = max(0, math_intercept)
        if x_intercept > np.min(x): x_intercept = np.min(x)

    return x_intercept, n_est

def robust_fit(circuit, freq, z, bounds, use_sigma=True):
    if circuit.initial_guess:
        lbs, ubs = bounds
        circuit.initial_guess = [max(lb+1e-12, min(g, ub-1e-12)) for g, lb, ub in zip(circuit.initial_guess, lbs, [u if u!=np.inf else 1e6 for u in ubs])]
    sigma = np.concatenate((np.abs(z), np.abs(z))) if use_sigma else None
    try:
        circuit.fit(freq, z, bounds=bounds, sigma=sigma)
        return True
    except:
        try:
            circuit.fit(freq, z, bounds=bounds, sigma=None)
            return True
        except: return False

def fit_semicircle_two_step(freq, z, min_idx):
    idx_cut = min_idx + 1
    f_semi, z_semi = freq[:idx_cut], z[:idx_cut]
    if len(z_semi) < 3: return None, []
    
    r0_est = np.min(np.real(z))
    if len(z_semi) > 0:
        r1_est = np.max(np.real(z_semi)) - r0_est
        idx_max_imag = np.argmax(-np.imag(z_semi))
        w_p = 2 * np.pi * f_semi[idx_max_imag]
    else:
        r1_est = np.max(np.real(z)) - r0_est
        w_p = 100
        
    c1 = CustomCircuit('R0-p(R1,CPE1)', initial_guess=[r0_est, r1_est, 1/(r1_est*w_p+1e-9), 0.9])
    if not robust_fit(c1, f_semi, z_semi, ([1e-6, 0, 1e-13, 0.5], [np.inf, np.inf, 1, 1.05]), True): return None, []
    
    r0, r1, q1, n1 = c1.parameters_
    _, n_tail = geometric_tail_analysis(z[min_idx:])
    
    # ==========================================================
    # 核心：寻优与误差计算解耦 (Decoupling Optimization from Error Estimation)
    # ==========================================================
    
    # 【准备第一跑】：构建“自适应补偿权重 (sigma_weighted)”，强行拉平采样密度分布，用于寻找最佳参数
    sigma_weighted = np.abs(z).copy()
    n_tail_pts = len(z) - min_idx
    if n_tail_pts > 0 and min_idx > 0:
        sigma_weighted[min_idx:] = sigma_weighted[min_idx:] / np.sqrt(max(1.0, (min_idx / n_tail_pts)))
        
    c_final = CustomCircuit('R0-p(R1,CPE1)-CPE2', initial_guess=[r0, r1, q1, n1, 1e-3, n_tail])
    
    try:
        # 【第一跑：寻优】：使用加权的 sigma_weighted 寻找全局最优参数 (防止高频密集点挟持)
        sigma_full_weighted = np.concatenate((sigma_weighted, sigma_weighted))
        c_final.fit(freq, z, bounds=([1e-6, 0, 0, 0.4, 1e-13, 0.15], [np.inf, np.inf, 1, 1.05, 1, 1.05]), sigma=sigma_full_weighted)
        
        # 【准备第二跑】：提取第一跑找到的最优参数作为新的初值
        best_params = c_final.parameters_
        
        # 重新定义电路（此时初值已经处于最优解坑底）
        c_stat_legal = CustomCircuit('R0-p(R1,CPE1)-CPE2', initial_guess=best_params)
        
        # 恢复原始的、未经篡改的模值权重，这是统计学合法性的根基
        sigma_original = np.abs(z).copy()
        sigma_full_original = np.concatenate((sigma_original, sigma_original))
        
        # 【第二跑：合法评估】：用真实的物理权重再拟合一次。
        # 因为初值已经完美，这步瞬间就能收敛。此次计算出的 conf_ 是 100% 统计学合法的。
        c_stat_legal.fit(freq, z, bounds=([1e-6, 0, 0, 0.4, 1e-13, 0.15], [np.inf, np.inf, 1, 1.05, 1, 1.05]), sigma=sigma_full_original)
        
        return c_stat_legal, ["R0", "R1", "Q1", "n1", "Q2", "n2"]
        
    except: 
        return None, []

def fit_tail_only(freq, z):
    rb_geom, n_geom = geometric_tail_analysis(z)
    hard_min_r0 = max(1e-9, np.min(np.real(z)) * 0.8)
    try:
        z_mod_low = np.abs(z[-1])
        f_low = freq[-1]
        q_est = 1.0 / (z_mod_low * (2 * np.pi * f_low)**n_geom + 1e-12)
        q_est = np.clip(q_est, 1e-13, 1.0)
    except:
        q_est = 1e-5 
        
    c = CustomCircuit('R0-CPE1', initial_guess=[max(rb_geom, hard_min_r0), q_est, n_geom])
    
    if not robust_fit(c, freq, z, ([hard_min_r0, 1e-13, 0.2], [np.inf, 1, 1.05]),  False): 
        return None, []
    return c, ["R0", "Q1", "n1"]

def trim_drift(f, z, min_idx):
    n = len(z)
    if n < 5 or min_idx >= n - 2:
        return f, z

    f_arc, z_arc = f[:min_idx+1], z[:min_idx+1]
    f_tail, z_tail = f[min_idx+1:], z[min_idx+1:]
    if len(z_tail) < 5: return f, z

    xr, yi = np.real(z_tail), -np.imag(z_tail)
    
    best_r2, best_slope, best_end_idx, best_start_idx = -1.0, 1.0, 0, 0
    window_size = max(5, int(len(xr) * 0.2))
    
    for i in range(0, len(xr) - window_size + 1):
        slope, _, r_val, _, _ = linregress(xr[i : i+window_size], yi[i : i+window_size])
        if r_val**2 > best_r2 and slope > 0.1:
            best_r2, best_slope, best_start_idx, best_end_idx = r_val**2, slope, i, i + window_size - 1

    ref_angle = np.rad2deg(np.arctan(best_slope))
    
    xr_gold, yi_gold = xr[best_start_idx : best_start_idx + window_size], yi[best_start_idx : best_start_idx + window_size]
    inner_dx, inner_dy = np.diff(xr_gold), np.diff(yi_gold)
    valid_mask = inner_dx > 0
    if np.any(valid_mask):
        inner_angles = np.rad2deg(np.arctan2(inner_dy[valid_mask], inner_dx[valid_mask]))
        sigma_theta = np.std(inner_angles)
        topology_weight = 1.5 if min_idx < 3 else 1.0 
        drift_tolerance = np.clip(topology_weight * (4.0 * sigma_theta + 5.0), 8.0, 20.0)
    else:
        drift_tolerance = 5.0 

    kept_indices = []
    xr_range = np.ptp(xr) if len(xr) > 0 else 1.0
    mono_threshold = xr_range * 0.005 
    for i in range(best_end_idx + 1):
        if i > 0 and xr[i] < xr[i-1] - mono_threshold: continue
        kept_indices.append(i)
        
    for i in range(best_end_idx + 1, len(xr)):
        dx, dy = xr[i] - xr[i-1], yi[i] - yi[i-1]
        if dx <= 0: continue 
        curr_angle = np.rad2deg(np.arctan2(dy, dx))
        if curr_angle >= ref_angle - drift_tolerance: kept_indices.append(i)
        else: break 

    if len(kept_indices) < len(xr) * 0.15: kept_indices = range(int(len(xr) * 0.4))
    keep_mask = np.array(kept_indices)
    return np.concatenate((f_arc, f_tail[keep_mask])), np.concatenate((z_arc, z_tail[keep_mask]))
      
    
def process_file(filepath, filename, plot_folder, output_pdf=False, img_dpi=300):
    print(f">>> 正在分析: {filename}", flush=True)
    f_raw, z_raw, f_rem, z_rem = load_and_clean_data(filepath)
    if f_raw is None or len(f_raw) < 10: return None
    
    z_backup_for_plot = z_raw.copy() 
    z_clean, bad_indices = filter_bad_z_prime(z_raw)
    if len(z_clean) < 10:
        print("   ❌ 物理清洗后数据点过少，文件作废")
        return None
        
    good_mask = np.ones(len(z_raw), dtype=bool)
    good_mask[bad_indices] = False
    f_raw = f_raw[good_mask]
    z_raw = z_raw[good_mask]
    
    # KK 校验 
    res_re, res_im, kk_ok = manual_kk_test(f_raw, z_raw)
    point_errs = np.maximum(np.abs(res_re), np.abs(res_im))
    mean_kk_error = np.mean(point_errs)
    
    # 形状初步预判
    global_shape, global_min_idx = detect_shape(z_raw)
    p_shape = "容抗弧+扩散尾" if global_shape == "semicircle_with_tail" else "单一直线特征"
    print(f"   🔎 初步识别形状: {p_shape} (Idx={global_min_idx})")
    fit_cache = {}

    # --核心逻辑区 ---
    def solve_impedance(gamma_val, current_g):
        valid_mask = point_errs <= gamma_val
        if np.sum(valid_mask) < 10: return None
        
        f_fit, z_fit = f_raw[valid_mask], z_raw[valid_mask]
        all_indices = np.arange(len(z_raw))
        kept_indices = all_indices[valid_mask]
        local_min_idx = np.argmin(np.abs(kept_indices - global_min_idx))
        
        try:
            drift_start = local_min_idx if global_shape == "semicircle_with_tail" else 0
            f_trimmed, z_trimmed = trim_drift(f_fit, z_fit, drift_start) 
            if len(f_trimmed) < 10: return None
            
            data_key = (len(f_trimmed), f_trimmed[0], f_trimmed[-1])
            if data_key in fit_cache: return fit_cache[data_key]

            # --- 步骤 A: 拟合基础模型 (直线) ---
            c_simple, p_simple = fit_tail_only(f_trimmed, z_trimmed)
            if not c_simple: return None
            z_pred_simple = c_simple.predict(f_trimmed)

            # --- 步骤 B: 探测候选点  ---
            scan_candidates = []
            track_msg = "" 

            if global_shape == "semicircle_with_tail":
                scan_candidates = [local_min_idx - 1, local_min_idx, local_min_idx + 1]
                track_msg = f"⚡ [快速微调] Idx={local_min_idx}"
            else:
                # ---单调性检查 (带平滑防抖) ---
                # 区分上升尾巴和半圆，且防止单点噪声误触发
                
                hf_limit = min(20, len(z_trimmed))
                if hf_limit > 5:
                    z_hf = z_trimmed[:hf_limit]
                    
                    # 1. 获取虚部数据
                    neg_imag = -np.imag(z_hf)
                    
                    # 平滑处理 (防抖)
                    # 窗口长度取 5 (或更小)，阶数取 1 (线性平滑) 或 2
                    # 抹平微小的噪声毛刺，使其回归单调
                    window_len = min(5, len(neg_imag))
                    if window_len % 2 == 0: window_len -= 1
                    if window_len >= 3:
                        try:
                            # polyorder=1 相当于滑动平均，适合判断单调性
                            neg_imag_smooth = savgol_filter(neg_imag, window_len, polyorder=1)
                        except:
                            neg_imag_smooth = neg_imag
                    else:
                        neg_imag_smooth = neg_imag

                    # 2. 寻找平滑后的峰值索引
                    peak_idx = np.argmax(neg_imag_smooth)
                    
                    # 3. 计算斜率翻转
                    xr_hf, yi_hf = np.real(z_hf), -np.imag(z_hf)
                    dx, dy = np.diff(xr_hf), np.diff(yi_hf)
                    slopes = dy / (dx + 1e-9)
                    sign_flips = np.sum(np.diff(np.sign(slopes)) != 0)
                    
                    # --- 核心裁决 ---
                    if sign_flips > 3:
                        track_msg = f"🚪 [拒绝:噪声] Flips={sign_flips}"
                    
                    # 【单调性判定】：检查平滑后的峰值位置
                    # 允许 1-2 个点的误差缓冲区 (如峰值在倒数第二个点，也算单调)
                    elif peak_idx >= len(neg_imag) - 2:
                        track_msg = f"🚪[拒绝:单调上升] 无半圆特征"
                    
                    else:
                        valley_candidate = np.argmin(slopes) + 1
                        scan_candidates = [valley_candidate, valley_candidate + 2] 
                        track_msg = f"🔍 [深度探查] 发现凸起 (Idx={peak_idx})"
                else:
                    track_msg = "🏃 [极速通道]"

            # 打印调试追踪 (仅在阈值变化较大时)
            if abs(gamma_val - current_g) < 0.001:
                print(f"      {track_msg}")

            # --- 步骤 C: 拟合复杂模型 (半圆) ---
            best_complex = None
            min_complex_rss = np.inf
            
            for s_idx in set(scan_candidates):
                if s_idx < 3 or s_idx >= len(z_trimmed) - 3: continue
                c_try, p_try = fit_semicircle_two_step(f_trimmed, z_trimmed, s_idx)
                if c_try and 0.2 <= c_try.parameters_[3] <= 1.1 and c_try.parameters_[1] > 0:
                    z_pred_try = c_try.predict(f_trimmed)
                    
                    # 寻找最佳候选时，同步加入 1.0 噪声底板
                    weight_try = np.abs(z_trimmed) + 1.0
                    rss_try = np.sum((np.abs(z_pred_try - z_trimmed) / weight_try)**2)
                    
                    if rss_try < min_complex_rss:
                        min_complex_rss = rss_try
                        best_complex = (c_try, p_try)

            # --- 步骤 D: F-Test 裁决 ---
            p_val = 1.0 
            final_c, final_p, current_shape = c_simple, p_simple, "tail_only"
            
            if best_complex:
                z_pred_complex = best_complex[0].predict(f_trimmed)
                p_val = calculate_f_test_pvalue(z_trimmed, z_pred_simple, z_pred_complex, 3, 6)
                
                # 高判定门槛 (99% 置信度)，防止高频伪影骗过 F-Test
                if p_val < 0.01:
                    final_c, final_p, current_shape = best_complex[0], best_complex[1], "semicircle_with_tail"

            # --- 步骤 E: 评分与结算 ---
            if current_shape == "semicircle_with_tail":
                rb = final_c.parameters_[0] + final_c.parameters_[1]
                if final_c.conf_ is not None and len(final_c.conf_) >= 2:
                    rb_std = np.sqrt(final_c.conf_[0]**2 + final_c.conf_[1]**2)
                else:
                    # [安全加固]：用预测残差反推置信度，加 1e-9 防止除零
                    z_pred_for_err = final_c.predict(f_trimmed)
                    rmse_rel = np.sqrt(np.mean((np.abs(z_trimmed - z_pred_for_err) / (np.abs(z_trimmed) + 1e-9))**2))
                    rb_std = rb * min(rmse_rel * 5.0, 1.0)
            else:
                rb = final_c.parameters_[0]
                if final_c.conf_ is not None and len(final_c.conf_) > 0:
                    rb_std = final_c.conf_[0]
                else:
                    z_pred_for_err = final_c.predict(f_trimmed)
                    rmse_rel = np.sqrt(np.mean((np.abs(z_trimmed - z_pred_for_err) / (np.abs(z_trimmed) + 1e-9))**2))
                    rb_std = rb * min(rmse_rel * 3.0, 1.0)

            u_rel = (rb_std / rb) * 100 if rb > 0 else 100
            retention_rate = np.sum(valid_mask) / len(z_raw)
            
           
            # 1. 线性惩罚 (1+gamma) 避免 丢失过多数据
            # 2. 留存率指数设为 3，鼓励保留更多数据
            current_score = (u_rel * (1.0 + gamma_val)) / (np.power(retention_rate, 3) + 1e-9)

            res = {
                "score": current_score, "circuit": final_c, "p_names": final_p, "rb": rb, 
                "r_err": u_rel, "shape": current_shape, "threshold": gamma_val,   
                "f_fit_final": f_trimmed, "z_final": z_trimmed, "valid_mask": valid_mask,
                "p_value": p_val
            }
            fit_cache[data_key] = res
            return res
        except: return None

    # --- 外层 Gamma 优化循环 ---
    median_e = np.median(point_errs)
    mad_e = np.median(np.abs(point_errs - median_e))
    # 上限到 5.0
    gamma_curr = np.clip(median_e + 2.5 * mad_e, 0.05, 5.0)
    
    step, max_iters, search_history, best_final = 0.08, 8, {}, None

    for _ in range(max_iters):
        candidates = [gamma_curr - step, gamma_curr, gamma_curr + step]
        iter_results = []
        for g in candidates:
           
            g = round(np.clip(g, 0.05, 5.0), 5)
            if g not in search_history:
                search_history[g] = solve_impedance(g, gamma_curr)
            if search_history[g]: iter_results.append(search_history[g])
        
        if not iter_results: break
        current_best = min(iter_results, key=lambda x: x['score'])
        if abs(current_best['threshold'] - gamma_curr) < 1e-6: step /= 2.5
        else: gamma_curr = current_best['threshold']
        best_final = current_best
        if step < 0.005: break 

    if best_final is None:
        print(f"   ❌ 算法收敛失败: {filename}"); return None

    # 打印最终决策信息
    f_pval = best_final["p_value"]
    f_shape = "容抗弧+扩散尾" if best_final['shape'] == "semicircle_with_tail" else "单一直线特征"
    print(f"   🏆 最佳阈值: {best_final['threshold']:.3f} | 最终形状: {f_shape} (P-val={f_pval:.2e})")
   
    # --- 统一解包并定义绘图变量 ---
    circuit = best_final["circuit"]
    rb = best_final["rb"]
    r_err = best_final["r_err"]
    shape = best_final["shape"]
    p_names = best_final["p_names"]
    valid_mask = best_final["valid_mask"]
    used_th = best_final["threshold"]
    
 
    f_fit_all = f_raw[valid_mask]   
    z_fit_all = z_raw[valid_mask]   
    f_fit_final = best_final["f_fit_final"] 
    z_final = best_final["z_final"]         
    
    if shape == "semicircle_with_tail":
        shape_str = "Arc+Tail"
        
        if circuit.conf_ is not None and len(circuit.conf_) >= 2: 
            rb_std = np.sqrt(circuit.conf_[0]**2 + circuit.conf_[1]**2)
        else: 
            rb_std = rb * 0.1
    else:
        shape_str = "TailOnly"
        
        if circuit.conf_ is not None and len(circuit.conf_) > 0: 
            rb_std = circuit.conf_[0]
        else: 
            rb_std = rb * 0.05
    rb_score = "Reliable" if r_err < 15 else "Review"

    # --- 绘图逻辑开始 ---
    fig = plt.figure(figsize=(20, 15)) 
    try:
      gs = GridSpec(3, 2, figure=fig, width_ratios=[1.3, 1], height_ratios=[1, 1, 1], hspace=0.25, wspace=0.15)
      ax1 = fig.add_subplot(gs[0:2, 0])
      ax2_fit = fig.add_subplot(gs[0, 1])
      ax2_kk = fig.add_subplot(gs[1, 1])
      ax2_zoom = fig.add_subplot(gs[2, 1])
      ax3 = fig.add_subplot(gs[2, 0])
      ax3.axis('off')

      if z_rem is not None and len(z_rem) > 0:
          ax1.plot(np.real(z_rem), -np.imag(z_rem), 'o', color='wheat', ms=6, label='Cleaned')
      if len(bad_indices) > 0:
          ax1.plot(np.real(z_backup_for_plot[bad_indices]), -np.imag(z_backup_for_plot[bad_indices]), 
                 'o', color='palegreen', ms=6, label='HF Pullback')
    
    # 绘制灰色 Trimmed 点
      truncated_mask = ~np.isin(f_fit_all, f_fit_final)
      if np.any(truncated_mask):
          ax1.plot(np.real(z_fit_all[truncated_mask]), -np.imag(z_fit_all[truncated_mask]), 'o', color='grey', ms=6, label='trimmed')
    
      if np.any(~valid_mask):
          ax1.plot(np.real(z_raw[~valid_mask]), -np.imag(z_raw[~valid_mask]), 'x', color='#696969', ms=5, alpha=0.9, label='KK Removed')
    
      ax1.plot(np.real(z_final), -np.imag(z_final), 'o', color='blue', ms=6, label='Fitted Data', zorder=5)
    
      f_line = np.logspace(np.log10(f_raw[-1]), np.log10(f_raw[0]), 200)
      z_pred = circuit.predict(f_line)
    
      ax1.plot(np.real(z_pred), -np.imag(z_pred), color='red', lw=2.5, label='Fit Model', zorder=6)
      ax1.plot(rb, 0, 'mX', ms=12, label=f'Rb={rb:.2f}Ω')
      ax1.axhline(0, color='k', lw=0.5)
      ax1.set_xlabel("Z' (Ohm)"); ax1.set_ylabel("-Z'' (Ohm)"); ax1.legend(loc='upper left'); ax1.grid(True, alpha=0.3); ax1.axis('equal')
      ax1.set_title(f"{filename}\nShape: {shape_str} | KK_Th: {used_th:.2f} | P-Val: {f_pval:.1e}", fontsize=13, fontweight='bold')

      z_p_final_pts = circuit.predict(f_raw)
      ax2_fit.plot(f_raw, (np.real(z_raw-z_p_final_pts)/np.abs(z_raw))*100, '.-')
      ax2_fit.axhline(0, color='black', lw=1, ls='--')
      ax2_fit.set_title("Fit Residuals (%)"); ax2_fit.set_xscale('log'); ax2_fit.grid(True, alpha=0.2)

      ax2_kk.plot(f_raw, res_re*100, 'c.-', label="KK Re%")
      ax2_kk.plot(f_raw, res_im*100, 'm.-', label="KK Im%")
      ax2_kk.axhline(0, color='black', lw=1, ls='--')
      ax2_kk.set_title(f"KK Consistency (Mean: {mean_kk_error*100:.2f}%)"); ax2_kk.set_xscale('log'); ax2_kk.grid(True, alpha=0.2)

      show_zoom = False
      if shape == "semicircle_with_tail":
          rs_val = circuit.parameters_[0]
          arc_dia = rb - rs_val
          range_x = np.max(np.real(z_raw)) - np.min(np.real(z_raw))
          range_y = np.max(-np.imag(z_raw)) - np.min(-np.imag(z_raw))
          max_span = max(range_x, range_y, 1e-6) 
          if arc_dia / max_span < 0.25:
              show_zoom = True
              ax2_zoom.plot(np.real(z_final), -np.imag(z_final), 'o', color='blue', ms=5)
              ax2_zoom.plot(np.real(z_pred), -np.imag(z_pred), color='red', lw=2)
              ax2_zoom.plot(rb, 0, 'mX', ms=10)
              pad = 0.25 * arc_dia
              ax2_zoom.set_ylim(-0.1 * arc_dia, 0.6 * arc_dia) 
              ax2_zoom.set_xlim(rs_val - pad, rs_val + arc_dia + pad)
              ax2_zoom.set_title("Semicircle Detail (Zoomed)", fontsize=11, color='darkblue')
              ax2_zoom.grid(True, alpha=0.2)
              ax2_zoom.set_aspect('equal', adjustable='datalim')
      if not show_zoom:
           ax2_zoom.text(0.5, 0.5, "No Detail View Required", ha='center', va='center', transform=ax2_zoom.transAxes, color='grey', fontsize=24, fontweight='bold')
           ax2_zoom.set_xticks([]); ax2_zoom.set_yticks([])

      def format_sci_latex(val):
          if 0.01 <= abs(val) <= 1000: return f"{val:.2f}"
          else:
              s = "{:.2e}".format(val)
              base, exp = s.split('e')
              return r"$" + base + r"\times 10^{" + str(int(exp)) + r"}$"

      t_vals = [[" Parameter", "Value ± 95% CI", "Error(%)", "Evaluation"]]
      kk_status = "Acceptable" if mean_kk_error < 0.08 else "Review"
      t_vals.append(["Overall KK Error", f"{mean_kk_error*100:.3f}%", "-", kk_status])
      rb_ci = rb_std * 1.96
      param_map = {"R0": "Rs (Ohm)", "R1": "Rion (Ion Migration)", "Q1": "CPE-Q", "n1": "CPE-n", "Q2": "Tail-Q", "n2": "Tail-n"}
      rb_row_idx = 2
      for i, n in enumerate(p_names):
          if shape == "tail_only" and n == "R0": continue 
          v = circuit.parameters_[i]
          s = circuit.conf_[i] if (circuit.conf_ is not None and len(circuit.conf_) > i) else 0
          display_name = param_map.get(n, n)
          per_err = (s/v)*100 if v!=0 else 0
          val_str = f"{format_sci_latex(v)} ± {format_sci_latex(s * 1.96)}"
          err_temp = format_sci_latex(per_err)
          if "$" in err_temp: err_str = err_temp[:-1] + r"\%$"
          else: err_str = f"{err_temp}%"
          eval_str = "Reliable" if per_err < 15 else "Review"
          t_vals.append([display_name, val_str, err_str, eval_str])

      tab = ax3.table(cellText=t_vals, loc='center left', bbox=[0, 0, 0.78, 1], colWidths=[0.25, 0.35, 0.25, 0.15])
      tab.auto_set_font_size(False); tab.set_fontsize(12)
      color_good, color_bad = '#C6EFCE', '#FFEB9C'
      for (row, col), cell in tab.get_celld().items():
          if row == 0: cell.set_facecolor('#D9D9D9'); cell.set_text_props(weight='bold')
          if col == 0: cell.set_text_props(ha='left') 
          else: cell.set_text_props(ha='center') 
          if col == 3 and row > 0:
              if row == rb_row_idx:
                  cell.set_facecolor(color_good if rb_score == "Reliable" else color_bad)
                  cell.get_text().set_text(" Pass" if rb_score == "Reliable" else " Review")
              else:
                  curr_eval = t_vals[row][3]
                  cell.set_facecolor(color_good if curr_eval in ["Reliable", "Excellent", "Acceptable"] else color_bad)
      tab.scale(1, 1.8)

      from matplotlib.patches import Rectangle
      rect = Rectangle((0.79, 0), 0.20, 1, transform=ax3.transAxes, linewidth=1, edgecolor='#D5D8DC', facecolor='#F8F9F9')
      ax3.add_patch(rect)
      criteria_text = ("$\mathbf{Scoring\ Criteria}$\n" " $\mathbf{Rb\ Error:}$\n" "   - Pass: < 15%\n" "   - Review: > 15%\n" " $\mathbf{Fit\ Error:}$\n" "   - Reliable: < 15%\n" "   - Review: > 15%\n" " $\mathbf{KK\ Consistency}$\n" "   - Acceptable: < 8%\n" "   - Review: > 8%")
      ax3.text(0.80, 0.5, criteria_text, transform=ax3.transAxes, va='center', ha='left', fontsize=14, linespacing=1.9)

      plt.subplots_adjust(left=0.05, right=0.97, top=0.94, bottom=0.04)
      if output_pdf:  # 如果界面勾选了输出PDF，才保存PDF
          plt.savefig(os.path.join(plot_folder, f"Result_{filename.split('.')[0]}.pdf"), bbox_inches='tight')
      plt.savefig(os.path.join(plot_folder, f"Result_{filename.split('.')[0]}.png"), dpi=img_dpi, bbox_inches='tight')
    finally:
      plt.close(fig)
    print(f"   ✅ 输出完成")
    return {"File": filename, "Rb": round(rb, 2), "Rb_Err%": round(r_err, 2), "Rb_Score": rb_score, "KK_Error%": round(mean_kk_error * 100, 3), "Used_KK_Th": round(used_th, 2), "Points_Rate": round(np.sum(valid_mask)/len(f_raw), 2)}

def run_main_program_entry(folder_path, output_pdf=False, img_dpi=300, progress_callback=None):
    print("="*50)
    print(f"🚀 单进程执行")
    print(f"📂 目标路径: {folder_path}")
    
    if not os.path.exists(folder_path):
        print(f"❌ 路径不存在: {folder_path}")
        return

    plot_folder = os.path.join(folder_path, "EIS分析图汇总")
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)

    all_files = [f for f in os.listdir(folder_path) if f.lower().endswith(('.txt', '.csv'))]
    if not all_files:
        print("❌ 错误: 没找到 .txt 数据文件！请检查文件格式。")
        return

    results = []
    for i, filename in enumerate(all_files):
        filepath = os.path.join(folder_path, filename)
        
        # 通知界面更新进度条
        if progress_callback:
            progress_callback(i, len(all_files), f"正在分析: {filename}")
            
        res = process_file(filepath, filename, plot_folder, output_pdf, img_dpi)
        if res:
            results.append(res)
            
    # 循环结束后，进度条拉满
    if progress_callback:
        progress_callback(len(all_files), len(all_files), "主程序分析完成！")

    if results:
        summary_path = os.path.join(folder_path, "EIS分析汇总表.csv")
        pd.DataFrame(results).to_csv(summary_path, index=False, encoding='utf-8-sig')
        print("="*50)
        print(f"✨ 处理完成！汇总表已生成。")
    else:
        print(f"❌ 路径不存在: {folder_path}")

# ==============================================================================
# 子程序 (Arrhenius Fit & Prediction) 
# ==============================================================================
import re
from itertools import combinations
kB = 8.617333262e-5
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Tahoma', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def parse_filename(filename):
    name_only = os.path.splitext(filename)[0]
    try:
        parts = re.findall(r'\d+', name_only)
        if len(parts) >= 3: return int(parts[0]), int(parts[1]), int(parts[2])
    except: return None, None, None
    return None, None, None



def get_best_valid_subsequence_sub(df_sample):
    df = df_sample.sort_values('Temp_C')
    n = len(df)
    if n < 3: return df.index.tolist() 
    
    # 内部辅助函数：计算最优折线（或直线）的残差平方和 RSS
    def get_min_rss(sub_df):
        X = sub_df['X'].values
        Y = sub_df['Y'].values
        m = len(X)
        if m <= 3:
            # 3个点以内只能拟合单一直线
            slope, intercept = np.polyfit(X, Y, 1)
            return np.sum((Y - (slope * X + intercept)) ** 2)
        
        best_rss = np.inf
        # 遍历所有可能的折点（保留两端至少包含2个点以构成线段）
        for i in range(1, m - 1):
            X_l, Y_l = X[:i+1], Y[:i+1]
            X_h, Y_h = X[i:], Y[i:]
            
            # 拟合低温段
            if len(X_l) >= 2: p_l = np.polyfit(X_l, Y_l, 1)
            else: p_l = [0, Y_l[0]]
            
            # 拟合高温段
            if len(X_h) >= 2: p_h = np.polyfit(X_h, Y_h, 1)
            else: p_h = [0, Y_h[0]]
            
            rss_l = np.sum((Y_l - (p_l[0]*X_l + p_l[1]))**2)
            rss_h = np.sum((Y_h - (p_h[0]*X_h + p_h[1]))**2)
            rss_total = rss_l + rss_h
            
            if rss_total < best_rss:
                best_rss = rss_total
        return best_rss

    best_indices = []
    max_len = 0
    best_rss_val = np.inf  # 残差
    
    # 限制遍历的最大长度，防止 combinations 导致的 O(2^n) 算力爆炸卡死
    search_limit = min(n, 18) 
    
    for r in range(search_limit, 2, -1): 
        for combo_indices in combinations(df.index, r):
            sub_df = df.loc[list(combo_indices)]
            
            # 物理约束1：必须单调递增（随温度升高，电导率必须上升）
            if not np.all(np.diff(sub_df['Y'].values) > 0): 
                continue
                
            # 物理约束2：整体趋势的斜率必须为负 (激活能 Ea > 0)
            slope, _, _, _, _ = linregress(sub_df['X'], sub_df['Y'])
            if slope >= 0: 
                continue 
            
            current_len = len(combo_indices)
            current_rss = get_min_rss(sub_df) # 计算该子集的最佳折线残差
            
            if current_len > max_len:
                max_len = current_len
                best_rss_val = current_rss
                best_indices = list(combo_indices)
            elif current_len == max_len:
                # 长度相同时，选残差(RSS)最小的
                if current_rss < best_rss_val:
                    best_rss_val = current_rss
                    best_indices = list(combo_indices)
                    
        # 只要找到了当前最长且符合物理规律的序列，就不再往下找更短的了
        if max_len > 0: 
            break
            
    return best_indices if best_indices else []
def fit_auto_breakpoint(df_sorted):
    """
    自动寻找最佳折线断点
    共享锚点(Shared Pivot) + 误差传递(Delta Method) + 点数安全检查
    """
    n_points = len(df_sorted)
    if n_points < 4: return None
    
    # 提取坐标
    X = df_sorted['X'].values
    Y = df_sorted['Y'].values
    
    # 获取 X 轴范围用于几何约束
    x_min, x_max = np.min(X), np.max(X)
    
    best_rss = np.inf
    best_model = None
    
    # 遍历每一个点作为“铰链点” (Pivot Point)
    # i 从 1 到 n-2，保证左右两边至少各有2个点
    for i in range(1, n_points - 1):
        
        # --- 1. 共享锚点切片 ---
        X_low = X[:i+1]  # 低温段包含锚点 i
        Y_low = Y[:i+1]
        X_high = X[i:]   # 高温段也包含锚点 i
        Y_high = Y[i:]
        
        # --- 2. 拟合与安全检查  ---
        # 只有当点数 > 2 时，才能计算协方差(cov=True)
        # 如果点数 == 2，误差无法定义(自由度为0)，只能算参数
        
        # 低温段拟合
        if len(X_low) > 2:
            try:
                (sl_l, int_l), cov_l = np.polyfit(X_low, Y_low, 1, cov=True)
            except: # 抓取矩阵奇异报错
                (sl_l, int_l) = np.polyfit(X_low, Y_low, 1)
                cov_l = np.zeros((2, 2))
        else:
            (sl_l, int_l) = np.polyfit(X_low, Y_low, 1)
            cov_l = np.zeros((2, 2))

        # 高温段拟合
        if len(X_high) > 2:
            try:
                (sl_h, int_h), cov_h = np.polyfit(X_high, Y_high, 1, cov=True)
            except:
                (sl_h, int_h) = np.polyfit(X_high, Y_high, 1)
                cov_h = np.zeros((2, 2))
        else:
            (sl_h, int_h) = np.polyfit(X_high, Y_high, 1)
            cov_h = np.zeros((2, 2))

        # --- 3. 物理与几何筛选 ---
        # 物理约束：斜率应为负 (Ea > 0)
        if sl_l >= 0 or sl_h >= 0: continue
        
        # 几何特征：低温斜率通常更陡 (更负)
        #if sl_l >= sl_h: continue (放宽限制适应更多的物质)

        # --- 4. 计算交点 ---
        denom = sl_l - sl_h
        if abs(denom) < 1e-9: continue
        x_int = (int_h - int_l) / denom
        
        # 几何约束：计算出的交点必须在当前锚点 i 附近
        # 允许一定容差 
        anchor_x = X[i]
        if abs(x_int - anchor_x) > (x_max - x_min) * 0.15:
            continue

        # --- 5. 计算残差平方和 (RSS) ---
        pred_l = sl_l * X_low + int_l
        pred_h = sl_h * X_high + int_h
        current_rss = np.sum((Y_low - pred_l)**2) + np.sum((Y_high - pred_h)**2)
        
        if current_rss < best_rss:
            best_rss = current_rss
            
            # --- 6. 误差传递  ---
            if len(X_low) > 2 and len(X_high) > 2:
                # 提取方差 (对角线)
                var_sl_l, var_int_l = cov_l[0,0], cov_l[1,1]
                var_sl_h, var_int_h = cov_h[0,0], cov_h[1,1]
                
                # 提取协方差 (非对角线：斜率与截距的负相关性)
                cov_m1_b1 = cov_l[0,1]
                cov_m2_b2 = cov_h[0,1]
                
                # 计算各自对交点 x 的方差贡献 
                var_x_l = var_int_l + (x_int**2) * var_sl_l + 2 * x_int * cov_m1_b1
                var_x_h = var_int_h + (x_int**2) * var_sl_h + 2 * x_int * cov_m2_b2
                
                # 总交点 x 的方差
                var_x = (var_x_l + var_x_h) / (denom**2)
                sigma_x = np.sqrt(max(0, var_x))
                
                # 换算成温度误差 dT = dx * 1000 / x^2
                tm_err = sigma_x * (1000.0 / (x_int**2))
            else:
                tm_err = 0.0 # 点数不足，无法估计误差

            calc_tm = 1000.0 / x_int - 273.15
            
             # 计算拟合系数 R2
            ss_tot_l = np.sum((Y_low - np.mean(Y_low))**2)
            r2_l = 1.0 - (np.sum((Y_low - pred_l)**2) / ss_tot_l) if ss_tot_l > 1e-12 else 1.0
            
            ss_tot_h = np.sum((Y_high - np.mean(Y_high))**2)
            r2_h = 1.0 - (np.sum((Y_high - pred_h)**2) / ss_tot_h) if ss_tot_h > 1e-12 else 1.0

            best_model = {
                'manual_tm': "Auto", 
                'calc_tm': calc_tm, 
                'tm_err': tm_err, 
                'intersection': (x_int, sl_l*x_int+int_l), 
                'low': (sl_l, int_l, r2_l, df_sorted.iloc[:i+1]), 
                'high': (sl_h, int_h, r2_h, df_sorted.iloc[i:])
            }

    # --- 如找不到折线，降级为单直线 ---
    if best_model is None:
        try:
            # 直线也计算协方差以保持格式一致
            if n_points > 2:
                (sl, int_c), cov = np.polyfit(X, Y, 1, cov=True)
            else:
                (sl, int_c) = np.polyfit(X, Y, 1)
                
            mid_x = (x_min + x_max) / 2.0
            calc_tm = 1000.0 / mid_x - 273.15
            
            best_model = {
                'manual_tm': "Auto (Linear)", 
                'calc_tm': calc_tm, 
                'tm_err': 0.0, # 单直线无折点误差
                'intersection': (mid_x, sl * mid_x + int_c), 
                'low': (sl, int_c, 0.99, df_sorted), 
                'high': (sl, int_c, 0.99, df_sorted)
            }
        except: pass
            
    return best_model

def fit_manual_breakpoint(df_sorted, manual_tm):
    """
    手动设定熔点断点
    使用强制约束的连续分段线性回归，确保两条线在设定的温度处完美交汇。
    """
    X = df_sorted['X'].values
    Y = df_sorted['Y'].values
    
    # 1. 将手动设定的摄氏度转为 1000/T 坐标 
    x_k = 1000.0 / (manual_tm + 273.15)
    
    # 2. 构建设计矩阵 (Design Matrix) 用于连续分段拟合
    # 数学方程: Y = y_k + m_low * X_low_shift + m_high * X_high_shift
    # 这确保了当 X = x_k 时，两边的 Y 值必定完全相等 (都等于 y_k)
    X_low_shift = np.where(X >= x_k, X - x_k, 0)
    X_high_shift = np.where(X < x_k, X - x_k, 0)
    
    A = np.vstack([np.ones_like(X), X_low_shift, X_high_shift]).T
    
    # 3. 最小二乘法求解最优参数 [y_k, m_low, m_high]
    # coef[0] = 交点处的 Y 值 (y_k)
    # coef[1] = 低温段斜率 (sl_l)
    # coef[2] = 高温段斜率 (sl_h)
    coef, _, _, _ = np.linalg.lstsq(A, Y, rcond=None)
    y_k, sl_l, sl_h = coef[0], coef[1], coef[2]
    
    # 4. 还原为普通直线方程的截距 
    # y = sl * x + int_c  =>  int_c = y_k - sl * x_k
    int_l = y_k - sl_l * x_k
    int_h = y_k - sl_h * x_k
    
    # 分割数据用于格式兼容
    mask_low = X >= x_k
    mask_high = X < x_k
    
    # 计算 R2
    pred = A.dot(coef)
    ss_tot = np.sum((Y - np.mean(Y))**2)
    r2 = 1 - np.sum((Y - pred)**2) / ss_tot if ss_tot > 0 else 0
    
    return {
        'manual_tm': manual_tm, 
        'calc_tm': manual_tm,     
        'tm_err': 0.0,              
        'intersection': (x_k, y_k), 
        'low': (sl_l, int_l, r2, df_sorted[mask_low]), 
        'high': (sl_h, int_h, r2, df_sorted[mask_high])
    }
   


def run_sub_program_entry(folder_path, thickness_um, diameter_mm, manual_tm_config, output_pdf=False, img_dpi=300, progress_callback=None):
    print(f"\n🔄 正在启动子程序 (Arrhenius Analysis)...")
    THICKNESS_CM = float(thickness_um) / 10000.0
    DIAMETER_CM = float(diameter_mm) / 10.0
    AREA_CM2 = np.pi * (DIAMETER_CM / 2)**2
    SIGMA_FACTOR = THICKNESS_CM / AREA_CM2 
    kB = 8.617333262e-5

    # .......................................................................................
    found_file = None
    for ext in[".xlsx", ".xls", ".csv"]:
        path = os.path.join(folder_path, f"EIS分析汇总表{ext}")
        if os.path.exists(path): found_file = path; break
    if not found_file: 
        print(f"❌ 错误：在路径 {folder_path} 下未找到 'EIS分析汇总表' 文件。")
        return
    print(f"Reading: {found_file}")
    df_summary = pd.read_csv(found_file) if found_file.endswith('.csv') else pd.read_excel(found_file)
    
    all_records = []
    for _, row in df_summary.iterrows():
        fname = str(row['File'])
        r, s, t = parse_filename(fname)
        status, reason = "Pending", ""
        sigma, x_val, y_val = 0, 0, 0
        if r is None: status, reason = "Rejected", "Filename Format Error"
        elif row['Rb'] <= 0: status, reason = "Rejected", "Rb <= 0"
        else:
            sigma = SIGMA_FACTOR / row['Rb']
            if not (1e-9 < sigma < 20.0): status, reason = "Rejected", "Conductivity Out of Range"
            else:
                x_val = 1000.0 / (t + 273.15)
                y_val = np.log10(sigma)
        all_records.append({'Ratio': r, 'Sample': s, 'Temp_C': t, 'File': fname, 'Rb': row['Rb'], 'Rb_Err': row.get('Rb_Err%', 5.0), 'Sigma': sigma, 'X': x_val, 'Y': y_val, 'Status': status, 'Reason': reason})
    
    full_df = pd.DataFrame(all_records)
    analysis_df = full_df[full_df['Status'] == 'Pending'].copy()
    output_dir = os.path.join(folder_path, "阿伦尼乌斯曲线图表结果")
    os.makedirs(output_dir, exist_ok=True)
    summary_report = []

    TARGET_TEMPS = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90]
    sigma_matrix_list = []      
    is_predicted_list = []      

    sorted_ratios = sorted(analysis_df['Ratio'].dropna().unique())
    for i, ratio in enumerate(sorted_ratios):
        if progress_callback:
            progress_callback(i, len(sorted_ratios), f"正在生成阿伦尼乌斯曲线: 序列 {ratio}")
            
        print(f"Analyzing Ratio: {ratio}wt%...")
        ratio_df = analysis_df[analysis_df['Ratio'] == ratio]
        
        # 使用界面传进来的 manual_tm_config，找不到默认用 Auto
        manual_tm = manual_tm_config.get(ratio, "Auto")
        
        valid_indices_all = []
        for sample_id in ratio_df['Sample'].unique():
            sample_df = ratio_df[ratio_df['Sample'] == sample_id]
            valid_idx = get_best_valid_subsequence_sub(sample_df)
            valid_indices_all.extend(valid_idx)
            rejected = list(set(sample_df.index) - set(valid_idx))
            full_df.loc[valid_idx, ['Status', 'Reason']] = ['Used (Valid)', 'Good']
            full_df.loc[rejected, ['Status', 'Reason']] = ['Rejected', 'Violates Trend']

        valid_df = full_df[(full_df['Ratio'] == ratio) & (full_df['Status'] == 'Used (Valid)')]
        invalid_df = full_df[(full_df['Ratio'] == ratio) & (full_df['Status'] == 'Rejected')]

        if valid_df.empty: 
            stats = pd.DataFrame(columns=['Temp_C', 'Y_mean', 'Y_std', 'X', 'Rb_Err'])
        else:
            stats = valid_df.groupby('Temp_C').agg({
                'Y': ['mean', 'std'], 'X': 'mean', 'Rb_Err': 'mean' 
            }).reset_index()
            stats.columns = ['Temp_C', 'Y_mean', 'Y_std', 'X', 'Rb_Err']
            stats['Y'] = stats['Y_mean']
            stats['Y_std'] = stats['Y_std'].fillna(0)

       
        # 双模并轨逻辑：计算值 vs 设定值
        auto_model = None
        current_model = None # 最终用于画图和预测的模型
        
        if not valid_df.empty:
            # 1. 第一跑：纯数学自动推断 
            auto_model = fit_auto_breakpoint(stats)
            
            # 2. 第二跑：根据配置决定画图模型
            if str(manual_tm).strip().lower() == "auto":
                current_model = auto_model
                manual_label_str = "Auto-Optimized"
            else:
                # 手动模式：强制在设定点交汇 
                current_model = fit_manual_breakpoint(stats, float(manual_tm))
                manual_label_str = f"{manual_tm} C"

        #将 current_model 赋值给 model
        model = current_model 

        row_data = {'Ratio': ratio}
        row_pred_status = {'Ratio': False}
        existing_temps_in_files = ratio_df['Temp_C'].unique()
        
        # 预测矩阵填充 
        for t_target in TARGET_TEMPS:
            match = stats[stats['Temp_C'] == t_target]
            final_sigma = None
            is_pred = False 
            if not match.empty:
                log_sigma = match['Y_mean'].values[0]
                final_sigma = 10 ** log_sigma
            elif t_target in existing_temps_in_files:
                if model:
                    x_target = 1000.0 / (t_target + 273.15)
                    sl_l, int_l, _, _ = model['low']
                    sl_h, int_h, _, _ = model['high']
                    # 分界点使用当前模型(Manual或Auto)的 calc_tm
                    boundary_tm = model['calc_tm'] 
                    pred_log = sl_l * x_target + int_l if t_target <= boundary_tm else sl_h * x_target + int_h
                    final_sigma = 10 ** pred_log
                    is_pred = True
            row_data[f"{t_target}C"] = final_sigma
            row_pred_status[f"{t_target}C"] = is_pred
            
        sigma_matrix_list.append(row_data)
        is_predicted_list.append(row_pred_status)

        # --- 绘图逻辑 ---
        fig, ax = plt.subplots(figsize=(10, 8))
        try:
            ax = plt.gca()
            if not invalid_df.empty: plt.scatter(invalid_df['X'], invalid_df['Y'], c='gray', s=30, alpha=0.4, label='Rejected')
            if not stats.empty: plt.errorbar(stats['X'], stats['Y_mean'], yerr=stats['Y_std'], fmt='o', color='blue', ecolor='blue', capsize=5, label='Valid Mean', zorder=5)

            res_txt = ""
            tm_err_str = ""
        
       
            if model:
                # 1. 取出用于画线的参数
                sl_l, int_l, r2_l, _ = model['low']
                sl_h, int_h, r2_h, _ = model['high']
                x_int, y_int = model['intersection']
            
            # 计算 Ea
                ea_l = -sl_l * 2.303 * kB * 1000
                ea_h = -sl_h * 2.303 * kB * 1000

            # 2. 取出用于显示的 Math Tm
                if auto_model:
                    math_tm_val = auto_model['calc_tm']
                    math_tm_err = auto_model['tm_err']
                else:
                    math_tm_val = model['calc_tm']
                    math_tm_err = model['tm_err']
            
                tm_err_str = f" ± {math_tm_err:.1f}" if math_tm_err < 10.0 else " (Uncertain)"

            # 3. 绘制辅助线
                x_split = x_int 
                if str(manual_tm).strip().lower() != "auto":
                    x_manual = 1000.0 / (float(manual_tm) + 273.15)
                    plt.axvline(x=x_manual, color='gray', linestyle='--', linewidth=2, label=f'Set Tm ({manual_tm}C)')
                    x_split = x_manual 

            # 4. 绘制拟合线
                x_start_l = stats['X'].max() if not stats.empty else x_split + 0.5
                plt.plot([x_start_l, x_split], [sl_l * x_start_l + int_l, sl_l * x_split + int_l], 'g-', lw=2.5, label='Low-T Fit')
             
                x_end_h = stats['X'].min() if not stats.empty else x_split - 0.5
                plt.plot([x_split, x_end_h], [sl_h * x_split + int_h, sl_h * x_end_h + int_h], 'orange', lw=2.5, label='High-T Fit')
            
                plt.scatter([x_int], [y_int], color='red', s=80, zorder=10, marker='D', label='Math Intersection')
            
             # 5. 生成对比文本
                res_txt = (f"Set Tm: {manual_label_str}\n"
                           f"Math Tm: {math_tm_val:.1f}{tm_err_str} C\n"
                           f"------------------\n"
                           f"Low-T Ea: {ea_l:.3f} eV\n"
                           f"High-T Ea: {ea_h:.3f} eV")
            
                summary_report.append({
                   'Ratio': ratio, 'Fit_Status': 'Success', 
                   'Manual_Set_Tm': manual_tm, 
                   'Math_Auto_Tm': round(math_tm_val, 2), 
                   'Tm_Error': round(math_tm_err, 2),
                   'Calc_Intersection_Tm_C': round(model['calc_tm'], 2),
                   'Ea_Low_eV': ea_l, 'R2_Low': r2_l, 
                   'Ea_High_eV': ea_h, 'R2_High': r2_h
               })
       
            
            # 绘制预测点
                valid_temps = stats['Temp_C'].values if not stats.empty else []
                missing_temps_to_plot = [t for t in existing_temps_in_files if t not in valid_temps]
                boundary_tm = model['calc_tm'] 
            
                for mt in missing_temps_to_plot:
                    mx = 1000.0 / (mt + 273.15)
                    my_pred = sl_l * mx + int_l if mt <= boundary_tm else sl_h * mx + int_h
                    plt.scatter([mx], [my_pred], c='gold', marker='X', s=150, edgecolors='k', zorder=10)
                    plt.text(mx, my_pred+0.05, f"{mt}C?", color='goldenrod', fontsize=9, ha='center')
            else:
                res_txt = "FIT FAILED\n(Slope > 0 or Not enough pts)"
                summary_report.append({'Ratio': ratio, 'Fit_Status': 'Failed'})

            
            plt.xlabel("1000 / T (1/K)")
            plt.ylabel("log10(Sigma) [S/cm]")
            plt.title(f"Arrhenius Plot: {ratio}wt% (Tm={manual_tm})", weight='bold')
            plt.grid(True, linestyle='--', alpha=0.5)
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            plt.legend(by_label.values(), by_label.keys(), loc='upper right', fontsize=9)
            props = dict(boxstyle='round', facecolor='white', alpha=0.9)
            plt.text(0.05, 0.05, res_txt, transform=ax.transAxes, verticalalignment='bottom', bbox=props, fontsize=10)
            if output_pdf:
                plt.savefig(os.path.join(output_dir, f"Arrhenius_{ratio}wt.pdf"), bbox_inches='tight')
            plt.savefig(os.path.join(output_dir, f"Arrhenius_{ratio}wt.png"), dpi=img_dpi, bbox_inches='tight')
        finally:    
            plt.close(fig) 


    excel_path = os.path.join(output_dir, "阿伦尼乌斯曲线分析表格信息.xlsx")
    if sigma_matrix_list:
        df_sigma_matrix = pd.DataFrame(sigma_matrix_list)
        df_is_predicted = pd.DataFrame(is_predicted_list)
        cols = ['Ratio'] + [c for c in df_sigma_matrix.columns if c != 'Ratio']
        df_sigma_matrix = df_sigma_matrix[cols]
        df_is_predicted = df_is_predicted[cols]
    else:
        df_sigma_matrix = pd.DataFrame()
        df_is_predicted = pd.DataFrame()

    try:
        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            pd.DataFrame(summary_report).to_excel(writer, sheet_name='Summary', index=False)
            if not df_sigma_matrix.empty:
                def highlight_predicted(data):
                    attr = 'background-color: #FFFF00' 
                    empty = ''
                    mask = df_is_predicted.values  
                    return pd.DataFrame(np.where(mask, attr, empty), index=data.index, columns=data.columns)
                (df_sigma_matrix.style
                 .apply(highlight_predicted, axis=None) 
                 .format("{:.2e}", subset=[c for c in cols if c != 'Ratio'], na_rep="")
                 .to_excel(writer, sheet_name='Sigma(Meas+Pred)', index=False))
            df_details = full_df[['Ratio', 'Sample', 'Temp_C', 'File', 'Rb', 'Sigma', 'Status', 'Reason']].copy()
            df_details.sort_values(by=['Ratio', 'Sample', 'Temp_C'], inplace=True)
            df_details.to_excel(writer, sheet_name='File Details', index=False)
        print(f"✅ 子程序分析完成！结果已保存至: {output_dir}")
    except Exception as e:
        print(f"❌ 保存Excel失败: {e}")
