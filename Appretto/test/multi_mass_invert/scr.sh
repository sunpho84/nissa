#!/bin/bash
(
    echo L 4
    echo T 8
    echo Seed 13
    echo TakeSlice 1
    echo TWall 0
    echo NoiseType 1
    echo Filename source
    echo Precision 64
) > input_generate_source

../../tools/generate_stochastic_source/generate_stochastic_source input_generate_source

(
    echo L 4
    echo T 8
    echo NMass 2
    echo Masses 0.50 0.60
    echo kappa 0.1770000
    echo GaugeConf ../Conf
    echo ThetaTXYZ 1 0 0 0
    echo Residue 1.e-6
    echo StoppingCriterion standard
    echo NiterMax 10000
    echo
    echo NSource 4
    echo
    echo Source source.00
    echo Output output.00
    echo
    echo Source source.01
    echo Output output.01
    echo 
    echo Source source.02
    echo Output output.02
    echo
    echo Source source.03
    echo Output output.03
) > input_multimass

../../tools/inverter/multimass_invert input_multimass

mkdir -p 0 1
for id in 0 1 2 3;do for im in 0 1; do for ud in 0 1;do mv output.0$id.0$im.$ud $ud/output.0$im.0$id; done; done; done

(
    echo L 4
    echo T 8
    echo TWall 0
    echo NPropFirstList 4
    echo 0/output.00 0.50 0 0 0
    echo 1/output.00 0.50 0 0 1
    echo 0/output.00 0.60 0 0 0
    echo 1/output.00 0.60 0 0 1
    echo NPropSecondList 4
    echo 0/output.00 0.50 0 0 0
    echo 1/output.00 0.50 0 0 1
    echo 0/output.00 0.60 0 0 0
    echo 1/output.00 0.60 0 0 1
    echo Ncontr 5
    echo 5 9
    echo 5 5
    echo 13 13
    echo 14 14
    echo 15 15
    echo NChromoContr 0
    echo Output two_points
) > input_contractions

../../projects/meson_2pts/meson_2pts input_contractions