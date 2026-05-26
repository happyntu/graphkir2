# ADR-004: C++23 改寫延後至 method 穩定後

- **狀態**：已接受
- **日期**：2026-05-26

## 背景

graphkir2 的計算瓶頸集中於兩個階段：

1. **EM allele typing** — 讀段與等位基因的似然矩陣反覆疊代，為 AVX2 SIMD 的理想目標
2. **CN depth accumulation** — BAM 逐位元計深度，適合 mmap 列式處理

`d:\works\bio\` 下的其他 C++23 專案（`haplophase` HMM forward-backward、`alignx` 列式 depth 查詢）已驗證相同模式的可行性。技術上可將熱路徑包成 `pybind11`/`nanobind` Python 擴充模組，外部工具呼叫（HISAT2、samtools、MUSCLE）與 pipeline 編排邏輯留在 Python。

## 決策

**延後 C++23 改寫，目前不執行。**

原因：

- 現階段首要目標是 3/5-digit 功能性精度（超越 Geny），不是執行速度
- synthetic stress sweep 仍有殘留問題（KIR2DS3 seed5102 一筆 5-digit 未解）
- 全面 C++23 改寫估計需 3–6 個月，會讓 method 工作完全停擺

## 後續行動

等以下條件成立後，再評估 C++23 extension：

1. method 面穩定（synthetic stress sweep 無殘留 candidate regression）
2. real-data HPRC 與 Geny 比較完成，paper-ready 結果確認
3. 以 **EM typing 似然矩陣計算** 為第一個 C++ extension 切入點，不動 guard 邏輯與 pipeline 結構

## 影響

- `src/graphkir2/` 繼續以純 Python 開發
- 效能優化暫以 `allele_base_top_n` gene-aware 設定（而非重寫）為主要手段
