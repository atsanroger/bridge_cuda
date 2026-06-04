# Float-Base Multiword-Precision Domain-Wall Solver
## 現況總結與「上超級電腦驗證」的提案

> 對象:教授面談（下週日）
> 目標:用目前在單張 RTX 3080 上得到的結果,論證下一步必須在超級電腦(A100/H100、多 GPU、production 格點)上測試。
> 日期:2026-06-04　平台:RTX 3080(SM 8.6, FP64 1:64, ~760 GB/s)、測試格點 32×8×8×12、Domain-Wall fermion、eo-preconditioned CGNR、mq=0.1。

---

## 0. 一句話(TL;DR)

我們在**只有 FP32 的消費級 GPU** 上,用 quasi multi-word 算術(Mukunoki–Ozaki, arXiv:2510.13536)實作出 Domain-Wall 解法器,**達到並超越原生 FP64 的精度**(triple-word 觀測量誤差 5.6×10⁻¹⁷ vs 原生 FP64 的 3.9×10⁻¹⁶)。效能上,所有保留意見都來自**小格點 + 單卡**這個測試平台的人為限制,而這些限制在 production 規模上**系統性地偏向對 multi-word 有利**。因此下一步應上超級電腦驗證。

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

**所有 §3 的效能保留意見,本質上都是「32×8×8×12 + 單卡」這個測試平台的人為限制,而且都系統性地偏向對 multi-word 不利。** production 規模會逆轉:

1. **快取假象 →(production)更純 memory-bound → MW overhead 下降。**
   3080 L2 只有 ~5MB,但 eo gauge ~10MB、測試 fermion 場很小 —— **測試格點有相當部分其實住在快取裡**。production 格點是 **GB 級**,什麼都進不了快取 → kernel 變得**更純粹 DRAM-bandwidth-bound** → MW 的三字元算術**更容易被延遲藏住** → 相對 overhead 應**往資料比 2×/3× 下降**(=論文 CPU 端「被 latency 藏住」的情形)。

2. **⚠️ §3.4 的 gauge 壓縮 −43% 很可能在大格點翻成正收益。**
   小格點上 gauge 讀取被快取吸收 → 省 bytes 省不到、白付算術 → −43%。大格點上 gauge 讀取是真的 DRAM 流量 → 省它划算 → **gauge-recon / 壓縮在 production size 可能變正的**。這條 lever**不該因小格點結果而放棄**。

3. **多 GPU → comm latency 主導 → MW overhead 相對更不顯著。**
   論文明言此點。我們 np=1 的數字是**最壞情況**;上規模後 MW 相對更便宜。

4. **更大體積 / 更小質量 → 條件數上升 → 迭代更多 → 捨入累積更多 → TW 精度領先擴大,且精度可能開始影響 nconv(論文的 fewer-iters 效應)。**

5. **A100/H100 → FP64 1:2 全速 → double-base QTW(~150-bit)回到論文 regime。**
   在 server 上做 **double-triple-word**,FP64 算術全速,於是「**modest overhead 換 beyond-FP64**」這個論文式賣點**會回來**。這是 3080 上做不到的(FP64 1:64)。

> 一句話:**我們量到的所有效能負面結論都偏保守/最壞情況;精度面的優勢則只會維持或擴大。** 要把這些「可能」變成數據,唯一的方法是上超算。

---

## 5. 提案:超級電腦上的實驗計畫（由便宜到昂貴)

| # | 實驗 | 平台 | 預測 / 要驗證的假說 | 成本 |
|---|---|---|---|---|
| 1 | **Volume scan**(遞增格點到超出快取) | 單 A100/H100 | MW overhead 隨體積**下降**趨向 2×/3×;recon −43% **翻正**;precision floor 隨體積走勢 | 低 |
| 2 | **Multi-GPU strong/weak scaling** | 多卡/多節點 | comm 主導下 MW overhead 相對**更小** | 中 |
| 3 | **double-base QTW** | A100/H100(FP64 1:2) | beyond-FP64 at modest overhead(論文 regime 重現) | 中 |
| 4 | **小夸克質量 / 物理 pion mass** | 同上 | 條件數爆炸區,FP64 也吃力;TW precision-reach + refinement 真正 engage | 高 |

**最小可行第一步 = 實驗 1(volume scan)**,單卡即可,直接證實/否證「小格點低估了我們」這個核心假說。

---

## 6. 我們在主張什麼 / 不主張什麼

- ✅ **主張**:在 FP32-only 硬體上達到並**超越 FP64 精度**(TW 5.6e-17 < FP64 3.9e-16),全程不碰 FP64。這是**無可取代**的能力(再多 FP64 也給不了 sub-double)。
- ✅ **主張**:DW 在沒有(好)FP64 的硬體上提供 **FP64 等級精度**;價值在 FP64 貧弱的 AI/消費級加速器。
- ✅ **主張**:此法是論文算法在 GPU/FP32/preconditioned-DWF 的首次落地,具發表價值。
- ❌ **不主張**(在 3080/單卡上):比 FP64 更快。論文的「fewer-iters 加速」對我們此設定不成立(nconv 相同),同基準 overhead 也輸論文。
- 🔬 **待驗證(=要上超算的理由)**:§4 的五點預測——production 規模 + 多 GPU + double-base 是否把效能劣勢逆轉成優勢。

---

## 附錄 A:可重現性

- 圖:`build/cuda-ampere-mpi/residual_compare.png`(由 `residual_floor_sweep.sh` 產生)。
- 腳本:`residual_floor_sweep.sh`(真實殘差掃描)、`compare_multiword.sh`(對照,已知遞迴殘差會誤導)、`ncu_profile.sh`(WSL2 不可用,留作 server 端用)。
- 量測方法注記:
  - 比較精度要用**真實殘差 ‖b−Ax‖**,不可用 CG 遞迴殘差(三精度重疊、且 1e-40 時皆 NaN breakdown)。
  - `max_restart=0` 會讓 CG 跑零迭代(m_Niter=Niter×Nrestart);需 ≥1。
  - bound 分析在 WSL2 用 `cuobjdump -res-usage` + `nvidia-smi dmon -s u`(ncu 詳細計數器在此 WSL2 不可用)。

## 附錄 B:關鍵數字速查

- Precision floor(真實殘差):FP32 1.6e-14 / DW ~1e-28 / TW ~3e-32。
- 觀測量精度:FP32 5e-10 / 原生 FP64 3.9e-16 / **TW 5.6e-17**。
- 牆鐘(vs FP64):DW ≈1.0× / TW ≈1.66×。同基準(vs FP32):DW 3.6× / TW 4.6×。
- TW hopb:167 reg + 2304B spill,~76% 峰值頻寬,dmon 90/68%。
- 已剔除 lever:recon(小格點 −43%,大格點待測)、shared-gauge(−14%)、occupancy(無益)、refinement@floor(零增益)。
