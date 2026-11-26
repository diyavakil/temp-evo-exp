1. Quality filtering: Preprocess raw ONT reads to get up to 100X cov using Filtlong
2. Align: Align filtered raw read to uploaded reference using minimap
3. Call variant: Call variants using Clair3 and bacteria tuned model
4. Filter variant: Filter low confidence and low depth variants
5. Annotate variant: Annotate and predict effect of variants using SnpEff
6. Subtract variants: Remove recurring artifacts between evolved lines at same temp
7. Filter candidates: keep protein-changing variants with VAF 0.9 or higher
