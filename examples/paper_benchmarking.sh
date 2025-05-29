# reset so no unintended runs
line_number=19
new_content="run_weather_loss = False"
sed -i "" "${line_number}s|.*|$new_content|" "bench/config_benchmark.ini"
line_number=21
new_content="random_seed = 24"
sed -i "" "${line_number}s|.*|$new_content|" "bench/config_benchmark.ini"
line_number=25
new_content="allocation_file = allocations/real_2018B.txt"
sed -i "" "${line_number}s|.*|$new_content|" "bench/config_benchmark.ini"
#
# # loop over the number of slots for each single shot request in the nominal request set with the real_2018B map
# outdir=bench/recreate_paper/ns_test/ns
# for ns_val in 12 6 3 2 1
# do
#   echo $outdir$ns_val/logfile-$ns_val.log
#   mkdir -p $outdir$ns_val/
#   echo "-ns " + $ns_val > $outdir$ns_val/logfile-$ns_val.log
#   astroq bench -nr 1 -ns $ns_val -cf bench/config_benchmark.ini > $outdir$ns_val/logfile-$ns_val.log
#   astroq plot -cf bench/config_benchmark.ini > $outdir$ns_val/logfile2-$ns_val.log
#   cp -R bench/outputs $outdir$ns_val/
#   cp -R bench/reports $outdir$ns_val/
# done
#
# # loop over the premade random allocation maps
# outdir=bench/recreate_paper/random_test
# for number in {0..9}
# do
#   echo $outdir/random_$number/logfile-$number.log
#   line_number=25
#   new_content="allocation_file = allocations/randoms/rand$number.txt"
#   sed -i "" "${line_number}s|.*|$new_content|" "bench/config_benchmark.ini"
#   echo "line number $line_number is now $new_content"
#   mkdir -p $outdir/random_$number/
#   astroq bench -nr 1 -ns 12 -cf bench/config_benchmark.ini > $outdir/random_$number/logfile-$number.log
#   astroq plot -cf bench/config_benchmark.ini > $outdir/random_$number/logfile2-$number.log
#   cp -R bench/outputs $outdir/random_$number/
#   cp -R bench/reports $outdir/random_$number/
# done
#
# # loop over the premade random every night single quarter allocation maps
# outdir=bench/recreate_paper/everynight_test
# for number in {0..9}
# do
#   echo $outdir/everynight_$number/logfile-$number.log
#   line_number=25
#   new_content="allocation_file = allocations/every_night/pattern_$number.txt"
#   sed -i "" "${line_number}s|.*|$new_content|" "bench/config_benchmark.ini"
#   echo "line number $line_number is now $new_content"
#   mkdir -p $outdir/everynight_$number/
#   astroq bench -nr 1 -ns 12 -cf bench/config_benchmark.ini > $outdir/everynight_$number/logfile-$number.log
#   astroq plot -cf bench/config_benchmark.ini > $outdir/everynight_$number/logfile2-$number.log
#   cp -R bench/outputs $outdir/everynight_$number/
#   cp -R bench/reports $outdir/everynight_$number/
# done

# reset so no unintended runs
line_number=25
new_content="allocation_file = allocations/real_2018B.txt"
sed -i "" "${line_number}s|.*|$new_content|" "bench/config_benchmark.ini"
# turn on weather loss
line_number=19
new_content="run_weather_loss = True"
sed -i "" "${line_number}s|.*|$new_content|" "bench/config_benchmark.ini"

# loop over the weather loss iterations with the real_2018B map
outdir=bench/recreate_paper/weather_test
starting_seed=24
for run_numb in {0..9}
do
  echo $outdir/weath_$run_numb/logfile-$run_numb.log
  mkdir -p $outdir/weath_$run_numb/
  line_number=21
  seed_value=$((run_numb + starting_seed))
  new_content="random_seed = $seed_value"
  sed -i "" "${line_number}s|.*|$new_content|" "bench/config_benchmark.ini"
  astroq bench -nr 1 -ns 12 -cf bench/config_benchmark.ini > $outdir/weath_$run_numb/logfile-$run_numb.log
  astroq plot -cf bench/config_benchmark.ini > $outdir/weath_$run_numb/logfile2-$run_numb.log
  cp -R bench/outputs $outdir/weath_$run_numb/
  cp -R bench/reports $outdir/weath_$run_numb/
done

# reset so no unintended runs
line_number=19
new_content="run_weather_loss = False"
sed -i "" "${line_number}s|.*|$new_content|" "bench/config_benchmark.ini"
line_number=21
new_content="random_seed = 24"
sed -i "" "${line_number}s|.*|$new_content|" "bench/config_benchmark.ini"
line_number=25
new_content="allocation_file = allocations/real_2018B.txt"
sed -i "" "${line_number}s|.*|$new_content|" "bench/config_benchmark.ini"
