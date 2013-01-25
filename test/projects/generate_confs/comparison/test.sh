#!/bin/bash

rm -fr conf output/*

$@ ../../../../projects/generate_confs/generate_stag ref_input
