rm -rf log.txt; touch log.txt

./run-all.bash f1 M1 z-positive >> log.txt &
./run-all.bash f1 M1 z-negative >> log.txt &
./run-all.bash f2 M1 z-positive >> log.txt &
./run-all.bash f2 M1 z-negative >> log.txt &
./run-all.bash f3 M1 z-positive >> log.txt &
./run-all.bash f3 M1 z-negative >> log.txt &
./run-all.bash f4 M1 z-positive >> log.txt &
./run-all.bash f4 M1 z-negative >> log.txt &
./run-all.bash f5 M1 z-positive >> log.txt &
./run-all.bash f5 M1 z-negative >> log.txt &
./run-all.bash f6 M1 z-positive >> log.txt &
./run-all.bash f6 M1 z-negative >> log.txt &
