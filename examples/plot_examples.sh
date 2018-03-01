#!/bin/bash
for s in `cat system_list`; do 
  nl=`cat ${s}.dat|wc -l`
  tail -n ${nl} ./example_${s}/Point.info|awk '{print $NF}' >./tmp.assign
  paste ${s}.dat ./tmp.assign >./tmp_assign.dat
  NCL=`grep TOTAL ./example_${s}/cluster.log |awk '{print $NF}'`
  rm ./tmp.assign
  echo "set term pngcairo size 600,1800" >tmp.gpl
  echo "set output 'plot_${s}.png'" >>tmp.gpl
  echo "set multiplot" >>tmp.gpl
  echo "set origin 0.0,0.66666" >>tmp.gpl
  echo "set size 1.0,0.33333" >>tmp.gpl
  if [ ${NCL} -eq 8 ];then 
    echo "set palette model RGB defined (0 'black',1 '#E9E202',2 '#D50AEC',3 '#F71513',4 '#1AF710',5 '#11E0F7',6 '#0216F8',7 'dark-green',8 '#F1D8EA')">>tmp.gpl
  fi
  if [ ${NCL} -eq 7 ];then 
    echo "set palette defined (0 'black',1 '#e41a1c',2 '#377eb8',3 '#4daf4a',4 '#984ea3',5 '#ff7f00',6 '#ffff33',7 '#a65628')">>tmp.gpl
  fi
  if [ ${NCL} -eq 3 ];then 
    echo "set palette defined ( 0 'black',1 'blue', 2 'yellow', 3 'red' )">>tmp.gpl
  fi
  if [ ${NCL} -eq 2 ];then 
    echo "set palette defined ( 0 'black',1 '#006400',  2 'blue' )">>tmp.gpl
  fi
  echo "unset colorbox;unset xtics;unset ytics;pl 'tmp_assign.dat' u 1:2:3 pt 7 ps 0.5 lc palette z t ''" >>tmp.gpl
  echo "reset" >>tmp.gpl
  echo "set origin 0.0,0.33333" >>tmp.gpl
  echo "set size 1.0,0.33333" >>tmp.gpl
  echo "unset colorbox;unset xtics;unset ytics" >>tmp.gpl
  if [ ${NCL} -eq 8 ];then 
    echo "set palette model RGB defined (1 '#E9E202',2 '#D50AEC',3 '#F71513',4 '#1AF710',5 '#11E0F7',6 '#0216F8',7 'dark-green',8 '#F1D8EA')">>tmp.gpl
  fi
  if [ ${NCL} -eq 7 ];then 
    echo "set palette defined (1 '#e41a1c',2 '#377eb8',3 '#4daf4a',4 '#984ea3',5 '#ff7f00',6 '#ffff33',7 '#a65628')">>tmp.gpl
  fi
  if [ ${NCL} -eq 3 ];then 
    echo "set palette defined ( 1 'blue', 2 'yellow', 3 'red' )">>tmp.gpl
  fi
  if [ ${NCL} -eq 2 ];then 
    echo "set palette defined ( 1 '#006400',  2 'blue' )">>tmp.gpl
  fi
  cat ./example_${s}/Dendrogram_graph_2.gpl >>tmp.gpl
  cp ./example_${s}/Dendrogram_labels_2.dat ./
  echo "reset" >>tmp.gpl
  echo "set origin 0.0,0.0" >>tmp.gpl
  echo "set size 1.0,0.33333" >>tmp.gpl
  cp ./example_${s}/Topography.mat ./
  cp ./example_${s}/mat_labels.dat ./
  echo "set palette grey negative;unset colorbox;unset border;unset xtics;unset ytics; plot 'Topography.mat' matrix with image t '','mat_labels.dat' u 1:2:3 w labels t '','' u 2:1:3 w labels t '' ">>tmp.gpl
  gnuplot tmp.gpl
  rm  ./tmp_assign.dat ./tmp.gpl ./Dendrogram_labels_2.dat ./Topography.mat ./mat_labels.dat
done
