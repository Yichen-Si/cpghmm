# HMM model for cumulated germline CpG methylation signature

### Inferred methylation level of autosomal CpGs
.bed files ready to view on UCSC Genome Browser is available [here](https://drive.google.com/drive/folders/11MQJvsCtgEOecii2u5X5gB-hkP2PIDVq?usp=share_link)

The scores are posterior probabilities of being hypo-methylated scaled to 0~1000. The estimates are based on TOPMed freeze 8 allele frequencies accessed through [Bravo](https://bravo.sph.umich.edu/freeze8/hg38/).


### Running on your own data

Dependence:
[Eigen](https://gitlab.com/libeigen/eigen)
[htslib](https://github.com/samtools/htslib)
[optim](https://github.com/kthohr/optim)

```sh
├── eigen
├── htslib
├── optim
│   └── build
│       ├── include
│       │   └── optim
│       └── lib
│           └── liboptim.so
└── cpghmm
```

(One might need to change the path to the dependencies in CMakeList.txt otherwise)

```sh
cd cpghmm
cmake .
make
```

Notes on installing optim:

```sh
# Move to the parent directory of cpghmm
git clone https://github.com/kthohr/optim ./optim
cd ./optim
export EIGEN_INCLUDE_PATH=path/to/eigen
mkdir -p build/include
mkdir build/lib
./configure -i "path/to/optim/build" -l eigen -p
make
make install
```

Example command
```sh
cpghmm cpg-cthmm --in-tsv ${input} --position_col ${x} --ac_col ${y} --sample-size ${n_sample} --init-emission ${emit} --init-transition ${tran} --out ${out} --region ${reg} --chunk-size ${ck_size} --optim-NM 1 --max-iter-outer 10 --max-iter-inner 100 --tol-outer 1e-5
# intput: tsv file containing a position column (specify by --position_col) and an allele count column (specify by --ac_col).
# emit, tran: (initial) emission and transition parameters. Example files are in /example
# out: prefix of output files
```

Output contains three files:

out.loo.tsv: leave-one-out likelihood

out.viterbi.tsv: viterbi path

out.likelihood: (marginal) conditional likelihood

If parameters are updated (--update-parameter 0 is not present), new parameter files are generated in the same format as the input initial parameter files.
