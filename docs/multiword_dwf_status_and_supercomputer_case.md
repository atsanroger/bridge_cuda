# Float-Base Multiword-Precision Domain-Wall Solver
## 現況總結與「上超級電腦驗證」的提案

> 對象:教授面談（下週日）
> 目標:用目前在單張 RTX 3080 上得到的結果,論證下一步必須在超級電腦(A100/H100、多 GPU、production 格點)上測試。
> 日期:2026-06-05　平台:RTX 3080(SM 8.6, FP64 1:64, ~760 GB/s)、Domain-Wall fermion、eo-preconditioned CGNR。
> 測試格點:32×8×8×12(mq=0.1)+ **新增 16⁴(Ns=8,三個夸克質量 0.10/0.05/0.025)** 過夜掃描。

---

## 0. 一句話(TL;DR)

我們在**只有 FP32 的消費級 GPU** 上,用 quasi multi-word 算術(Mukunoki–Ozaki, arXiv:2510.13536)實作出 Domain-Wall 解法器,**達到並超越原生 FP64 的精度**(triple-word 觀測量誤差 5.6×10⁻¹⁷ vs 原生 FP64 的 3.9×10⁻¹⁶),且此精度優勢在 16⁴ / 三質量上**不隨體積改變**(已驗)。
效能面**誠實說**:我們的瓶頸是 hopping(Deo)kernel 的 **register pressure**,而它**與體積無關、也不會被多卡的 comm 重疊藏掉**(我們的設計 comm 與 Deo 重疊 → 全程 Deo-bound)。因此原先「大格點/多卡會自動把 overhead 壓下來」的樂觀預期**已被自己的資料修正掉**。效能面真正站得住的論點只剩**一個**:在 **FP64 全速的 A100/H100 上做 double-base QTW**,回到論文的 modest-overhead regime。**這就是要上超算的核心理由**(外加在 production 體積驗證 precision-reach)。
**而且這個主打實驗已經不是「提案」—— 2026-06-05 已在 3080 上實證會動、會贏**:double-base double-triple 多字在真正的 FP64 算術上,真實殘差地板達 ~8.6×10⁻⁷⁶,**比原生 FP64 低約 46 個數量級**,且明確勝過 double-double ~13 個數量級(見 §2.4)。3080 只是太慢(FP64 1:64),要的就是 A100/H100 把它跑快。

---

## 1. 我們建了什麼

- 在 Bridge++ 的 alt_Accel(CUDA)後端,為 eo-preconditioned Domain-Wall 5din 解法器加入三條精度路徑,以 YAML 切換(`multiword_mode: FP | DW | TW`):
  - **FP** = 單精度 float(baseline)
  - **DW** = quasi double-word(float-float,≈2⁻⁴⁸,~雙精度等級)
  - **TW** = quasi triple-word(float-triple,≈2⁻⁷²,**超越** FP64)
- 全程 **FP32-only**:device 端從不使用 FP64(常數在 host 端用 double-double 預備,device 端維持 float-float/float-triple)。
- 涵蓋:hopping 算符、5D 預條件子(LU/Dinv)、CG 的 BLAS 層(norm2/dot/axpy/aypx,各自的 dual-/triple-word reduction)、以及全 forward 算符 M 與其 adjoint(供 refinement 使用)。

> 學術定位:這是論文 quasi-MW 算法**第一個** GPU / FP32 / eo-preconditioned-DWF 的落地,直接回答該論文 Future Work 點名的三個 open question(preconditioning、FP32-base、GPU)。

---

## 2. 主結果:Precision Reach —— 在 FP32-only 硬體上贏過 FP64

### 2.0 一張圖看完整條精度階梯(附圖 `sweep16_6mode_20260605_2156/ladder_mq0.10.png`,三個質量各一張)

**16⁴(Ns=8)、三個夸克質量、同一套 eo-preconditioned CGNR**,只用 YAML 的 `field_type` + `multiword_mode` 切換,六條精度同台比較。真實殘差地板(mq=0.10):

| base | 模式 | 真實殘差地板 | 備註 |
|---|---|---|---|
| FP32(消費級 GPU) | single | 3.3×10⁻¹³ | float 上限 |
| FP32 | float-float (QDW) | 9.5×10⁻²⁸ | ~FP64 級 |
| FP32 | **float-triple (QTW)** | **6.9×10⁻³⁴** | **← 在 FP32 硬體上就贏過原生 FP64** |
| FP64(伺服器 GPU) | native | 3.1×10⁻³⁰ | FP64 上限 |
| FP64 | double-double (QDW) | 2.2×10⁻⁶² | |
| FP64 | **double-triple (QTW)** | **3.2×10⁻⁷⁶** | |

兩個關鍵交叉都在圖上可見:**(a)** FP32-base 的 float-triple 跌破原生 FP64(消費級 GPU 主張);**(b)** FP64-base 的 double-triple 到 ~10⁻⁷⁶(伺服器主張,§2.4)。

**三質量(0.10 / 0.05 / 0.025)趨勢 ── 兩個小質量的重要發現**(`ladder_mq0.05.png`、`ladder_mq0.025.png`):
- 質量越小、條件數越大 → CG 越慢 → 要更多迭代才觸底。給足迭代(cap 2400)後 double-triple 在 mq=0.025 反而到更深的 **~1.9×10⁻⁹⁰**(double-double ~3.3×10⁻⁶²、native FP64 ~4×10⁻³⁰)。
- **FP32 base 在小質量觸不到地板**:float run 設 cap=550(刻意停在 FP32 遞迴殘差崩潰點 ~600 之前)。mq=0.025 時 float DW/TW 跑到 cap 仍在 ~6×10⁻¹⁹ **下降中、尚未觸底**,而它再多跑就會崩 → 等於「float 在小質量根本到不了自己的精度地板」;反觀 FP64-base multiword 一路下探到 1e-90。這是「小夸克質量該用 double-base / 伺服器」的具體論據。(圖中 FP32 虛線到 cap 後以點線水平延伸代表「保持最後值」,非新迭代。)

> 註:這些殘差圖只解 12 個 color-spin source 的**第一個**(殘差下降曲線只需一個),故**不產生有效的 2pt 物理觀測量**;觀測量請看 §2.2 的完整收斂跑。另有一張單格點 4⁴×8 的同款圖 `multiword_precision_ladder.png` 作為交叉驗證。

下面 §2.1–2.4 拆開細談。



### 2.1 真實殘差 ‖b−Ax‖ 隨迭代下降(附圖 `residual_compare.png`)

掃描迭代上限,記錄**真實殘差**(非 CG 遞迴殘差;後者三精度完全重疊、不反映精度)。三精度在 ~150 迭代前完全一致(**CG 收斂率與精度無關**),之後各自觸底:

| 迭代上限 | FP32 | DW (QDW) | TW (QTW) |
|---|---|---|---|
| 150 | 2.3e-13 | 2.2e-13 | 2.2e-13 |
| 200 | **1.3e-14** ←觸底 | 5.1e-16 | 5.1e-16 |
| 325 | 1.58e-14 | 7.8e-21 | 8.1e-21 |
| 500 | 1.58e-14 | 1.4e-28 | **3.1e-32** |

- **FP32 硬觸底於 ‖b−Ax‖ ≈ 1.6×10⁻¹⁴**,且低於原生 FP64 的結果精度。
- DW / TW 往下鑽 **14+ 個數量級**;TW 在迭代後段(>400)再拉開 DW 約 4 個數量級。

### 2.2 收斂後的物理觀測量精度(criterion 1e-30,完整 2pt)

| 模式 | nconv | 觀測量誤差 | 說明 |
|---|---|---|---|
| FP32 | 495 | ~5×10⁻¹⁰ | 單精度上限 |
| 原生 FP64 | (參考) | ~3.9×10⁻¹⁶ | FP64 硬體上限 |
| DW | 495 | ~FP64 等級 | dual-word floor |
| **TW** | 495 | **5.6×10⁻¹⁷** | **比 FP64 準 ~7×** |

> 關鍵:**精度不改變收斂次數**(三者皆 nconv=495),multi-word 的價值是把**可達精度地板**壓低,不是減少迭代。

### 2.3 16⁴ 三質量過夜掃描(2026-06-05,附圖 `overnight_16p4_20260605_0821/residual_mq*.png`)

把格點放大到 **16⁴(Ns=8)**、掃三個夸克質量,用內建 `true_residual_interval` 一次求解即得完整 true-residual 曲線。**精度結論與小格點完全一致,且不隨體積改變**:

| mq | nconv(FP/DW/TW) | FP32 floor | DW floor | TW floor | 觀測量:DW vs TW | FP 誤差 |
|---|---|---|---|---|---|---|
| 0.10  | 484 / 484 / 484 | 3.3e-13 | 9.5e-28 | **1.5e-30** | 一致到 ~3e-16(FP64 級) | ~3e-7 |
| 0.05  | 741 / 741 / 741 | 5.9e-13 | 1.5e-27 | **1.0e-30** | 一致到 ~3e-16 | ~3e-7 |
| 0.025 | 875 / 874 / 874 | 6.1e-13 | 1.4e-27 | **1.8e-30** | 一致到 ~3e-16 | ~3e-7 |

- **三精度迭代數相同**(每質量),且隨質量變小而增加(484→741→875,條件數上升)——再次證實**精度買的是地板、不是迭代**。
- TW 在迭代後段(>~780)再次拉開 DW 約 2–3 個數量級(eps³ 在 DW 逼近 float-float floor 後才顯威)。`residual_mq0.025.png` 是最乾淨的代表圖。
- **物理觀測量:DW 與 TW 互相一致到 ~3×10⁻¹⁶(FP64 級),FP32 偏離 ~3×10⁻⁷** —— 與小格點同樣的故事,在 production-relevant 體積上重現。

### 2.4 ★ DOUBLE-base 多字:在 FP64 硬體上超越原生 FP64(2026-06-05 已實證)

**這是給超級電腦提案的關鍵新證據。** 把 base 型別從 float 換成 **double**(`field_type: double_eo`),同一套 quasi-MW 演算法就變成 **double-double(QDW)/ double-triple(QTW)**,在**真正的 FP64 算術**上跑。

> **16⁴ 已驗證(§2.0 的 `ladder_mq*.png`)**:double-triple 在 mq=0.10 到 3.2×10⁻⁷⁶、mq=0.025 到 **1.9×10⁻⁹⁰**(double-double ~3×10⁻⁶²、native FP64 ~3×10⁻³⁰)── 即在 production-relevant 體積上,FP64-base 多字穩定超越原生 FP64 **~45–60 個數量級**。下表為單格點 4⁴×8 的精準驗證(附圖 `double_base_residual_compare.png`):

| 模式 | 真實殘差地板 ‖b−Ax‖²/‖b‖² | 相對原生 FP64 |
|---|---|---|
| FP(原生 FP64) | 1.59×10⁻³⁰ | —(這就是 FP64 的極限) |
| **DW(double-double)** | **1.03×10⁻⁶²** | **低 ~32 個數量級** |
| **TW(double-triple)** | **8.57×10⁻⁷⁶**(且在迭代上限仍在下降) | **低 ~46 個數量級** |

- **在 FP64 硬體上,double-triple 多字比原生 FP64 多約 46 個數量級;且 TW 明確勝過 DW ~13 個數量級**(三字真的比兩字準,不是被夾住)。三者 nconv=202 全收斂、零錯誤、FP 路徑仍精確命中參考值 4.35052428395743e-01。
- 這正是論文 regime 的 GPU 版:**modest 額外成本換 beyond-FP64**,且映射到 A100/H100(FP64 1:2 全速)。**等於把「上超算」的主打實驗,先在 3080 上證明會動、會贏** —— 只是 3080 的 FP64 1:64 讓它慢(故 timing 無意義,要的是 A100)。
- 實作:(1) 修好 `ConstantMemoryTraits<double>` 的 mid/lo stub + 補上 `const_a_double`(原本為 0,導致 double+TW 卡在 1e-11 不收斂);(2) 加入 host **triple-double(~159-bit)** 常數管線(`compute_precond_coeffs_td`,建在既有 double-double 原語上),把 double-base QTW 從「被夾在 double-double」提升到 genuine double-triple。改動 2 檔:`constant_memory.cu`、`constant_memory_inline.h`(+ `validate_double_tw.sh`)。FP native-double 路徑不受影響。

### 2.5 ★ 真實(非 tiled)config 驗證 + 我們和 Ozaki 的 regime 差異(2026-06-06)

**為什麼要做這節**:先前所有測試的 gauge 是把一個 4⁴×8 block 週期性自我複製的 tiled config(見 `src/lib/IO/fieldIO_Text_4x4x4x8.cpp`),其 Dirac 算符可 Bloch 分解 → 條件數人工地良態。最大的疑慮是:**「精度只影響 floor、不影響迭代」會不會只是 tiling 假象?換成真 config 是否會出現 Ozaki 那種「低精度多跑迭代」效應?** 我們用 Bridge 內建 quenched HMC(β=6.0)在 CPU 上自產了一個**真實 thermalized 16⁴ config**(plaq≈0.596,equilibrium),重跑三質量 × 六精度。

**結論:descent overlap 是真的物理,不是 tiling 假象。** 附圖 `sweep16_real_q450/ladder_real_mq*.png`:

| mode | mq=0.10 floor | mq=0.05 | mq=0.025 | 到 1e-12 的迭代數 |
|---|---|---|---|---|
| float FP32 | 3.0e-13 | 3.7e-13 | 3.7e-13 | **140**(所有精度相同)|
| float-float QDW | 7.5e-28 | 8.0e-28 | 1.0e-27 | 140 |
| float-triple QTW | 5.1e-37 | 2.3e-31 | 1.3e-28 | 140 |
| double FP64 | 2.4e-30 | 3.4e-30 | 3.4e-30 | 140 |
| double-double QDW | 1.8e-62 | 1.9e-62 | 2.9e-62 | 140 |
| **double-triple QTW** | **6.5e-90** | **1.7e-90** | **1.1e-90** | 140 |

- 真 config 上,**六種精度全部在 iter 140 抵達 1e-12**,完全不按精度分離(tiled 同質量是 160,真 config 其實**略快**)。floor 階梯與 tiled 幾乎一致 → **精度決定 floor、spectrum 決定 rate 的分工,在真實 gauge 上依然成立**。
- tiled caveat 解除一半:tiling 確實讓 nconv / volume-scaling 不 physical(那部分仍須謹慎),但 **descent overlap 本身是真現象**,不是自我複製造成的。

**我們和 Mukunoki–Ozaki(2510.13536)為何差那麼多 —— 兩個便宜的驗證**

差異的唯一旋鈕是**條件數 κ**,不是精度本身。機制是 finite-precision CG 的 orthogonality loss:其發作早晚 ∝ ε × 迭代數 × κ。我們靠 **eo-preconditioning + 良態 DWF** 把自己鎖在小 κ,所以多精度買到的是 **reach(floor)** 而非 **speed(迭代)**;Ozaki 在無預條件的病態 SuiteSparse SPD(κ up to 1e10)上,FP32 早早失去正交性 → 多精度才會**減少迭代**。兩個獨立證據:

1. **驗證 1 —— 真 Bridge log 的 recursive vs true residual**(附圖 `verify1_recursive_vs_true.png`)。log 的 `[TRUERES]` 記了兩欄:recursive(CG 內部遞迴)與 true ‖b−Ax‖。六種精度的 **recursive descent 完全疊合(跨精度離散度全程 < 0.18 dex,不到 1.5 倍)**,只有 **true residual 的 floor 隨 ε 階梯排開**。這就是「rate 與精度無關、precision 只設 reach」的硬證據。

2. **驗證 3 —— 脫離 Bridge 的 toy CG κ-sweep**(附圖 `toy_cg_kappa_sweep.png`,純 numpy、無預條件):

   | κ | FP64 到 1e-10 | FP32 到 1e-10 | FP32 floor |
   |---|---|---|---|
   | 10 | 18 iter | **18 iter(完全相同)** | 9e-11 |
   | 10³ | 166 iter | **永不到(stall)** | 4e-10 |
   | 10⁵ | 1197 iter | **永不到** | 7e-6 |
   | ≥10⁷ | 兩者皆失敗 | — | ~1 |

   κ=10(我們的 regime)FP32/FP64 迭代數**完全相同**;κ≥10³(Ozaki regime)FP32 因失去正交性**收斂不了**。同一段程式碼,精度對 rate 的影響**完全由 κ 開關** → 證實機制理解無誤。

**對提案的意涵(兩面誠實講)**:✅ 我們的 **reach 主張很硬**且已在真 config 證實(floor 1e-90,比 FP64 深 ~60 個數量級,Ozaki 做不到此深度);⚠️ 我們**不主張**省迭代 —— 在 eo-preconditioned DWF 上多精度不省迭代,已用真 config 直接證實(140 iter,六精度全同)。要復現 Ozaki 的「省迭代」效應,須走進病態 regime(關/弱化 preconditioner、或 mq→物理點產生近零模)── 這正是下一個值得在超算上做的實驗(見 §5)。

---

## 3. 效能:單張 RTX 3080 上的誠實評估

### 3.1 牆鐘時間(代表值,完整 2pt)

| 對照軸 | DW | TW |
|---|---|---|
| **vs 原生 FP64**(產品目標軸) | ≈ **1.0×**(FP64 等級精度、FP64 等級速度、**不需 FP64 硬體**) | ≈ **1.66×**(但**精度超越 FP64**) |
| vs FP32 同基準(算法效率軸) | 3.6× | 4.6× |
| 論文同基準(FP64-base, CPU) | QDW 1.3× | QTW 2.4× |

**誠實面對**:以「MW vs 自己的 base」這把尺,我們(FP32-base)比論文(FP64-base)差。原因見 §3.3。

### 3.2 kernel 分解(TW)

- **4D hopping stencil(`bulk`,等同論文 SpMV)佔 TW 求解 ~74%** —— 是主成本。
- 5D 預條件子(`Dee_inv`)僅 2–3× vs FP32(已優化到接近 memory-bound)。
- CG 的 BLAS/normalization 僅約 ~12%。

### 3.3 Bound 分析(無 ncu;WSL2 不支援 ncu 詳細計數器,改用 `cuobjdump` + `nvidia-smi dmon`)

| | FP32 hopb | TW hopb |
|---|---|---|
| Registers / thread | 72 | **167** |
| Local-memory spill | 0 | **2304 B** |
| dmon sm% / mem%(整個 solve) | 69 / 51 | 90 / 68 |
| 達峰值 DRAM 頻寬 | — | ~76% |

**結論:bandwidth-bound,但被 register pressure 壓住 occupancy**(167 reg + 2304B spill → 低 occupancy → 記憶體延遲藏不滿 → 只到 ~76% 峰值頻寬)。**不是 ALU/compute-bound**。

### 3.4 已實測、可剔除的優化迷思(避免教授問起)

- **SU(3) 重建(gauge 18→12)在小格點 = −43% 慢**:重建的三字元算術沒被藏住。**但這很可能是小格點假象**(見 §4)。
- **Shared-memory gauge staging = −14%**(register-bound kernel 不吃這套)。
- **提高 occupancy(大 block / color-split)在 float 上無增益**:已貼著頻寬天花板。
- **Iterative refinement 在 criterion 1e-30 = 零增益**(TW 已在 floor;只有**小夸克質量**才 engage)。

---

## 4. 為什麼單卡小格點的效能數字「低估」了這個方法（提案核心）

**所有 §3 的效能保留意見,本質上都是「小格點 + 單卡」這個測試平台的人為限制。** 但 16⁴ 掃描帶來一個**誠實的修正**(見下方 ⚠️),反而讓「必須上多卡 / double-base」的結論更硬:

> **⚠️ 16⁴ 實測修正(2026-06-05):單卡放大體積『不足以』降低 overhead。**
> 32×8×8×12 → 16⁴(~2GB)後,solver-only overhead **沒有**下降,維持 **DW ~3.8× / TW ~5.9× vs FP32**(三質量皆然,per-iter 成本固定)。原因:16⁴ 在 3080 上 **仍是 register/ILP-bound**(bulk float-triple stencil 佔 TW 求解 63–66%,167 reg + spill 壓住 occupancy),**還沒**進入純 DRAM-bound regime。**結論:下方第 1 點的「體積→overhead 下降」crossover,單卡 16⁴ 達不到 —— 需要(a)遠更大的體積、(b)多 GPU(comm 主導)、或(c)double-base(FP64 全速)。這正是要上超算的數據佐證,而非反證:單卡單軸放大已被證實不夠。**

1. **快取假象 →(production)更純 memory-bound → MW overhead 下降。** —— 方向仍成立,但 **16⁴ 還不夠大**(見上方 ⚠️):需到真正塞滿 DRAM 流量的體積才會 crossover。
   3080 L2 只有 ~5MB,但 eo gauge ~10MB、測試 fermion 場很小 —— **小測試格點有相當部分其實住在快取裡**。production 格點是 **GB 級**,什麼都進不了快取 → kernel 變得**更純粹 DRAM-bandwidth-bound** → MW 的三字元算術**更容易被延遲藏住** → 相對 overhead 應**往資料比 2×/3× 下降**(=論文 CPU 端「被 latency 藏住」的情形)。16⁴(~2GB)是這條曲線上的一個下界資料點,證明 crossover 在更大體積處。

2. **⚠️ §3.4 的 gauge 壓縮 −43% 很可能在大格點翻成正收益。**
   小格點上 gauge 讀取被快取吸收 → 省 bytes 省不到、白付算術 → −43%。大格點上 gauge 讀取是真的 DRAM 流量 → 省它划算 → **gauge-recon / 壓縮在 production size 可能變正的**。這條 lever**不該因小格點結果而放棄**。

3. **多 GPU → comm 主導 → MW 更便宜?── 此論點對我們的設計『不成立』(已自我修正)。**
   我們的 hopping 設計**把 halo 通訊與 bulk Deo 計算重疊**(comm overlapped with Deo)。只要 Deo-compute > comm-time,通訊就被完全藏住、kernel 維持 **Deo-bound**。MW 反而**放大** Deo,使 comm 更容易被藏 → runtime 全程跟著 Deo 走 → **MW/FP overhead 比值被『保持』,不會因多卡而縮小**。先前實驗也證實多 GPU **仍是 Deo 主導、非 comm 主導**。
   ⚠️ 另外:**QDW/QTW 目前『尚未』實作多 GPU 支援** —— 多卡是一個**開發前置工作**,不是現成可跑的實驗。所以多卡掃描的價值是「驗證 overhead 比值是否如預期被保持/或在 FP 變 comm-bound 的極高節點數時才略縮」,而**不是**原先寫的「comm 主導讓 MW 變便宜」。

4. **更大體積 / 更小質量 → 條件數上升 → 迭代更多 → 捨入累積更多 → TW 精度領先擴大,且精度可能開始影響 nconv(論文的 fewer-iters 效應)。**

5. **A100/H100 → FP64 1:2 全速 → double-base QTW(~150-bit)回到論文 regime。**
   在 server 上做 **double-triple-word**,FP64 算術全速,於是「**modest overhead 換 beyond-FP64**」這個論文式賣點**會回來**。這是 3080 上做不到的(FP64 1:64)。

> 一句話:**我們量到的所有效能負面結論都偏保守/最壞情況;精度面的優勢則只會維持或擴大。** 要把這些「可能」變成數據,唯一的方法是上超算。

---

## 5. 提案:超級電腦上的實驗計畫（由便宜到昂貴)

| # | 實驗 | 平台 | 預測 / 要驗證的假說 | 成本 |
|---|---|---|---|---|
| 1 | **Volume scan**(遞增格點到**遠超** 16⁴/快取) | 單 A100/H100 | 16⁴ 單卡已證 overhead **未**降(仍 register-bound);要找出 overhead 真正趨向 2×/3× 的 crossover 體積;recon −43% **翻正**點 | 低 |
| ★ | **double-base QTW(主打)** | A100/H100(FP64 1:2) | FP64 全速 → beyond-FP64 at modest overhead(論文 regime 重現)。**這是效能面唯一站得住的論點** | 中 |
| 2 | **病態 regime:小夸克質量→物理 pion mass / 弱化 preconditioner** | A100/H100 | **驗證 2(復現 Ozaki)**:升高 κ → 測 recursive descent 是否按精度「裂開」(低精度多跑迭代)。§2.5 已用 toy CG 證實機制由 κ 開關、用真 16⁴ config 證實當前良態 regime 下不裂;此實驗找出 DWF 的裂開臨界,讓「省迭代」主張可成立或明確否決 | 高 |
| 3 | (前置)**多 GPU 支援開發 + scaling** | 多卡/多節點 | 先**實作** QDW/QTW 的多卡(目前無);驗證 overhead 比值在 Deo-overlap-comm 設計下**被保持**(非 comm 主導) | 高(含開發) |

**最小可行第一步 = ★ double-base QTW**(單張 A100/H100 即可),直接證實「FP64-fast 硬體上 MW 回到論文的 modest-overhead regime」。
volume scan(實驗 1)仍值得做,但 16⁴ 單卡已顯示體積放大**未必**降 overhead(register-bound),故期望要放低。多 GPU 是**開發前置**、不是現成實驗。

---

## 6. 我們在主張什麼 / 不主張什麼

- ✅ **主張**:在 FP32-only 硬體上達到並**超越 FP64 精度**(TW 5.6e-17 < FP64 3.9e-16),全程不碰 FP64。這是**無可取代**的能力(再多 FP64 也給不了 sub-double)。
- ✅ **主張**:DW 在沒有(好)FP64 的硬體上提供 **FP64 等級精度**;價值在 FP64 貧弱的 AI/消費級加速器。
- ✅ **主張**:此法是論文算法在 GPU/FP32/preconditioned-DWF 的首次落地,具發表價值。
- ❌ **不主張**:比 FP64 更快。論文的「fewer-iters 加速」對我們此設定不成立(nconv 相同),同基準 overhead 也輸論文。
- ❌ **不再主張(已自我修正)**:「大格點 / 多卡會自動降低 overhead」。16⁴ 證實單卡放大體積**未**降 overhead(register-bound、與體積無關);且我們 comm 與 Deo 重疊 → 全程 **Deo-bound**,多卡**不會**讓 comm 主導來相對壓低 MW。**這兩條預測作廢。**
- 🔬 **待驗證(=要上超算的『主要』理由)**:**double-base QTW 在 FP64 全速硬體(A100/H100)上重現論文的 modest-overhead beyond-FP64 regime**。次要:在 production 體積上驗證 precision-reach 維持。
- 🛠️ **前置工程**:QDW/QTW **尚無多 GPU 支援**,多卡實驗需先開發。

---

## 附錄 A:可重現性

- 圖:`build/cuda-ampere-mpi/residual_compare.png`(由 `residual_plot.sh` 產生)。
- **原生功能(已內建)**:YAML `Solver: true_residual_interval: N` —— CG 解法器每 N
  迭代記錄**真實殘差** ‖b−Ax‖(tagged `[TRUERES]`,在欄位自身的 multiword 精度下重算,
  FP32-only)。所以**單次求解**即得完整 true-residual 下降曲線,不必再用掃描法。
- 腳本:`residual_plot.sh`(原生 `[TRUERES]`,**每模式只跑一次**,產生上圖)、
  `residual_floor_sweep.sh`(舊版迭代上限掃描,24 次,留作對照)、
  `compare_multiword.sh`(已知遞迴殘差會誤導,留作反例)、
  `ncu_profile.sh`(WSL2 不可用,留作 server 端用)。
- 量測方法注記:
  - 比較精度要用**真實殘差 ‖b−Ax‖**,不可用 CG 遞迴殘差(三精度重疊、且 1e-40 時皆 NaN breakdown)。
  - `max_restart=0` 會讓 CG 跑零迭代(m_Niter=Niter×Nrestart);需 ≥1。
  - bound 分析在 WSL2 用 `cuobjdump -res-usage` + `nvidia-smi dmon -s u`(ncu 詳細計數器在此 WSL2 不可用)。

## 附錄 B:關鍵數字速查

- Precision floor(真實殘差):FP32 ~1.6e-14 / DW ~1e-28 / TW ~3e-32(32×8×8×12);16⁴ 同序:FP32 ~3–6e-13 / DW ~1e-27 / TW ~1–1.8e-30。
- 觀測量精度:FP32 5e-10 / 原生 FP64 3.9e-16 / **TW 5.6e-17**。16⁴:DW≈TW 一致到 ~3e-16,FP32 偏 ~3e-7。
- 牆鐘(vs FP64):DW ≈1.0× / TW ≈1.66×。同基準(vs FP32):32³格點 DW 3.6×/TW 4.6×;**16⁴ solver-only DW ~3.8×/TW ~5.9×(未隨體積下降——仍 register-bound)**。
- nconv 與精度無關:32×8×8×12 三者=495;16⁴ 三者=484/741/875(隨 mq↓ 上升)。
- TW hopb:167 reg + 2304B spill,~76% 峰值頻寬,dmon 90/68%;16⁴ bulk 佔 TW 求解 63–66%。
- 已剔除 lever:recon(小格點 −43%,大格點待測)、shared-gauge(−14%)、occupancy(無益)、refinement@floor(零增益)。
- 16⁴ 過夜掃描產物:`docs/overnight_16p4_20260605_0821/`(3 張 residual PNG、summary.tsv、MASTER.log、9 份 run log)。
