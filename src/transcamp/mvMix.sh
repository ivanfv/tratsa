#/bin/bash

benchs=("result_power-MPIII-SVF_n180000_w1325_s1_d0_t256_he8_hm23_le8_lm7_*.csv"
        "result_power-MPIII-SVF_n180000_w1325_s1_d0_t256_he8_hm23_le5_lm10_*.csv"
        "result_power-MPIII-SVF_n180000_w1325_s1_d0_t256_he8_hm23_le5_lm2_*.csv"
        "result_power-MPIII-SVF_n180000_w1325_s1_d0_t256_he8_hm7_le5_lm10_*.csv"
        "result_power-MPIII-SVF_n180000_w1325_s1_d0_t256_he8_hm7_le5_lm2_*.csv"
        "result_audio-MPIII-SVD_w200_s1_d0_t256_he8_hm23_le8_lm7_*.csv"
        "result_audio-MPIII-SVD_w200_s1_d0_t256_he8_hm23_le5_lm10_*.csv"
        "result_audio-MPIII-SVD_w200_s1_d0_t256_he8_hm23_le5_lm2_*.csv"
        "result_audio-MPIII-SVD_w200_s1_d0_t256_he8_hm7_le5_lm10_*.csv"
        "result_audio-MPIII-SVD_w200_s1_d0_t256_he8_hm7_le5_lm2_*.csv"
        "result_seismology-MPIII-SVE_n180000_w50_s1_d0_t256_he8_hm23_le8_lm7_*.csv"
        "result_seismology-MPIII-SVE_n180000_w50_s1_d0_t256_he8_hm23_le5_lm10_*.csv"
        "result_seismology-MPIII-SVE_n180000_w50_s1_d0_t256_he8_hm23_le5_lm2_*.csv"
        "result_seismology-MPIII-SVE_n180000_w50_s1_d0_t256_he8_hm7_le5_lm10_*.csv"
        "result_seismology-MPIII-SVE_n180000_w50_s1_d0_t256_he8_hm7_le5_lm2_*.csv"
        "result_e0103_n180000_w500_s1_d0_t256_he8_hm23_le8_lm7_*.csv"
        "result_e0103_n180000_w500_s1_d0_t256_he8_hm23_le5_lm10_*.csv"
        "result_e0103_n180000_w500_s1_d0_t256_he8_hm23_le5_lm2_*.csv"
        "result_e0103_n180000_w500_s1_d0_t256_he8_hm7_le5_lm10_*.csv"
        "result_e0103_n180000_w500_s1_d0_t256_he8_hm7_le5_lm2_*.csv"
        "result_penguin_sample_TutorialMPweb_w800_s1_d0_t256_he8_hm23_le8_lm7_*.csv"
        "result_penguin_sample_TutorialMPweb_w800_s1_d0_t256_he8_hm23_le5_lm10_*.csv"
        "result_penguin_sample_TutorialMPweb_w800_s1_d0_t256_he8_hm23_le5_lm2_*.csv"
        "result_penguin_sample_TutorialMPweb_w800_s1_d0_t256_he8_hm7_le5_lm10_*.csv"
        "result_penguin_sample_TutorialMPweb_w800_s1_d0_t256_he8_hm7_le5_lm2_*.csv"
        "result_human_activity-MPIII-SVC_w120_s1_d0_t256_he8_hm23_le8_lm7_*.csv"
        "result_human_activity-MPIII-SVC_w120_s1_d0_t256_he8_hm23_le5_lm10_*.csv"
        "result_human_activity-MPIII-SVC_w120_s1_d0_t256_he8_hm23_le5_lm2_*.csv"
        "result_human_activity-MPIII-SVC_w120_s1_d0_t256_he8_hm7_le5_lm10_*.csv"
        "result_human_activity-MPIII-SVC_w120_s1_d0_t256_he8_hm7_le5_lm2_*.csv")

destDir="results/results_mixed_allHigh_profileLow"
mkdir $destDir

for i in "${benchs[@]}"; do
  echo "$i --> $destDir"
  mv results/$i $destDir
done;

