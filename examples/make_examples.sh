for s in `cat system_list`; do 
  mkdir example_${s}
  cp ${s}.dat ./example_${s}
  cp repeat_${s}.in ./example_${s}
  cd ./example_${s}
  ../../bin/DPA.x -ignze <repeat_${s}.in>${s}.std.out
  cd ..
done
