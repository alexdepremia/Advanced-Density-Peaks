for s in `cat system_list`; do 
  mkdir example_${s}_scan
  cp ${s}.dat ./example_${s}_scan
  cp ${s}_scan.in ./example_${s}_scan
  cd ./example_${s}_scan
  ../../bin/DPA_scan.x -ignze <${s}_scan.in>${s}.std.out
  cd ..
done
