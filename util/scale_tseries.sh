sed 's/^/10*/' epg.txt | bc | awk '{printf "%f\n", $0}' > epg_scaled_10.txt

