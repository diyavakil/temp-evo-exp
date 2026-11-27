All code below optimized to be run in Google Colab.

Estimate ANI using skani:
```
# install skani
%%bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
./bin/micromamba create -y -n skani_env -c conda-forge -c bioconda skani
```

```
# run skani comparison between two genomes
%%bash
./bin/micromamba run -n skani_env skani dist \
  "[GENOME_1].fasta" \ # make sure [GENOME_1].fasta is uploaded to the runtime
  "[GENOME_2].fasta" # make sure [GENOME_1].fasta is uploaded to the runtime
```


1. Quality filtering: Preprocess raw ONT reads to get up to 100X cov using Filtlong
2. Align: Align filtered raw read to uploaded reference using minimap
3. Call variant: Call variants using Clair3 and bacteria tuned model
4. Filter variant: Filter low confidence and low depth variants
5. Annotate variant: Annotate and predict effect of variants using SnpEff
6. Subtract variants: Remove recurring artifacts between evolved lines at same temp
7. Filter candidates: keep protein-changing variants with VAF 0.9 or higher
