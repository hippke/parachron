S=all_mn1_g3_dt180_nsub2_un_ga120
for ep in "$@"; do s=$((ep-100)); [ $s -lt 0 ] && s=0; python3 slice_err.py $S $s 100 cal_$ep.npy; done
