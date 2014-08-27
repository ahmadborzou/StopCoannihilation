#!/bin/bash

hadd results_PhaseII4_B_14TEV_HT1_140PileUp.root `ls results_PhaseII4_B_14TEV_HT1_140PileUp_*`;
for i in {1..5}; do hadd results_PhaseII4_BB_14TEV_HT${i}_140PileUp.root `ls results_PhaseII4_BB_14TEV_HT${i}_140PileUp_*`; done;
for i in {1..3}; do hadd results_PhaseII4_BBB_14TEV_HT${i}_140PileUp.root `ls results_PhaseII4_BBB_14TEV_HT${i}_140PileUp_*`; done;
for i in {1..4}; do hadd results_PhaseII4_H_14TEV_HT${i}_140PileUp.root `ls results_PhaseII4_H_14TEV_HT${i}_140PileUp_*`; done;
for i in {1..6}; do hadd results_PhaseII4_LL_14TEV_HT${i}_140PileUp.root `ls results_PhaseII4_LL_14TEV_HT${i}_140PileUp_*`; done;
for i in {1..3}; do hadd results_PhaseII4_LLB_14TEV_HT${i}_140PileUp.root `ls results_PhaseII4_LLB_14TEV_HT${i}_140PileUp_*`; done;
for i in {1..5}; do hadd results_PhaseII4_TB_14TEV_HT${i}_140PileUp.root `ls results_PhaseII4_TB_14TEV_HT${i}_140PileUp_*`; done;
for i in {1..5}; do hadd results_PhaseII4_TJ_14TEV_HT${i}_140PileUp.root `ls results_PhaseII4_TJ_14TEV_HT${i}_140PileUp_*`; done;
for i in {1..4}; do hadd results_PhaseII4_TTB_14TEV_HT${i}_140PileUp.root `ls results_PhaseII4_TTB_14TEV_HT${i}_140PileUp_*`; done;
for i in {1..4}; do hadd results_PhaseII4_BJJ_14TEV_HT${i}_140PileUp.root `ls results_PhaseII4_BJJ_14TEV_HT${i}_140PileUp_*`; done;








