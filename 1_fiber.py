import matplotlib.pyplot as plt
from matplotlib import font_manager
import math

# === 字體設定（避免中文亂碼） ===
plt.rcParams['font.sans-serif'] = ['Microsoft JhengHei']
plt.rcParams['axes.unicode_minus'] = False

print("=== 多模光纖脈衝展寬模擬器 ===")

# === 基本輸入 ===
n1 = float(input("請輸入纖芯折射率 n1 (例如 1.480): "))
n2 = float(input("請輸入包層折射率 n2 (例如 1.460): "))
c = 3e8  # 真空光速 (m/s)
L_start = float(input("請輸入起始長度 (公里): "))
L_end = float(input("請輸入結束長度 (公里): "))
L_step = float(input("請輸入間隔 (公里): "))

# === 光纖與光源參數 ===
a = float(input("請輸入纖芯半徑 a (μm，例如 25): ")) * 1e-6  # μm → m
lam = float(input("請輸入工作波長 λ (nm，例如 1300): ")) * 1e-9  # nm → m
tau0 = float(input("請輸入初始脈衝寬度 τ0 (ns，例如 0.1): ")) * 1e-9  # ns → s

# === 進階輸入（材料色散與光源寬度） ===
delta_lambda = input("輸入光源波長展寬 Δλ (nm) (預設 0.2): ")
dispersion_coeff = input("輸入材料色散係數 D (ps/(nm·km)) (預設 17): ")

delta_lambda = float(delta_lambda) if delta_lambda else 0.0
dispersion_coeff = float(dispersion_coeff) if dispersion_coeff else 0.0

# === 類型選擇 ===
fiber_type = input("光纖型態：1=Step-Index, 2=Graded-Index (預設 1): ")
fiber_type = fiber_type.strip()
fiber_type = 2 if fiber_type == '2' else 1

# === 計算初步參數 ===
Delta = (n1**2 - n2**2) / (2 * n1**2)       # 相對折射率差
NA = math.sqrt(n1**2 - n2**2)                # 數值孔徑
V = (2 * math.pi * a / lam) * NA             # 歸一化頻率
M = V**2 / 2                                 # 模式數近似

print("\n=== 光纖參數 ===")
print(f"相對折射率差 Δ = {Delta:.6f}")
print(f"數值孔徑 NA = {NA:.4f}")
print(f"歸一化頻率 V = {V:.3f}")
print(f"估計模式數 M ≈ {M:.1f}")

# === 儲存模擬結果 ===
L_list = []
delta_t_total_list = []
tau_out_list = []

# === 模擬 ===
L = L_start
while L <= L_end:
    L_m = L * 1000  # km → m

    # (1) 模式色散 Δt_modal
    if fiber_type == 1:  # Step-index
        delta_t_modal = (n1 * Delta * L_m) / c
    else:  # Graded-index
        delta_t_modal = (n1 * (Delta ** 2) * L_m) / (8 * c)

    # (2) 材料色散 Δt_material
    delta_t_material = (dispersion_coeff * delta_lambda * L) * 1e-12  # ps→s

    # (3) 波導色散 Δt_waveguide（若已包含於 D 可設為 0）
    delta_t_waveguide = 0

    # (4) 合併（RMS）
    delta_t_rms = math.sqrt(delta_t_modal**2 +
                            delta_t_material**2 +
                            delta_t_waveguide**2)

    # (5) 最終脈衝寬度（含輸入脈衝 τ0）
    tau_out = math.sqrt(tau0**2 + delta_t_rms**2)

    # 儲存結果（轉成 ns）
    L_list.append(L)
    delta_t_total_list.append(delta_t_rms * 1e9)
    tau_out_list.append(tau_out * 1e9)

    L += L_step

# === 輸出結果表格 ===
print("\n=== 模擬結果 ===")
print(" 長度 (km) │ 模式+色散展寬 Δt_total (ns) │ 輸出脈衝寬度 τ_out (ns)")
print("---------------------------------------------------------------")
for L, dt, tout in zip(L_list, delta_t_total_list, tau_out_list):
    print(f"{L:>10.1f} │ {dt:>27.6f} │ {tout:>21.6f}")

# === 視覺化 ===
plt.figure(figsize=(9, 5))
plt.plot(L_list, delta_t_total_list, 'o-', label="Δt_total (展寬)", color='royalblue')
plt.plot(L_list, tau_out_list, 's--', label="τ_out (輸出脈衝寬度)", color='darkorange')
plt.title("多模光纖脈衝展寬模擬 (Pulse Broadening vs Fiber Length)")
plt.xlabel("光纖長度 L (km)")
plt.ylabel("脈衝展寬 (ns)")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()
