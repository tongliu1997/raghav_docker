source /home/tl543/raghav_docker/setup.sh
input_list="testfile.list"
out_name="test.root"
system="Ru"
if [ $# -ge 1 ]; then
  input_list=$1
fi
if [ $# -ge 2 ]; then
  out_name=$2
fi
if [ $# -ge 3 ]; then
  system=$3
fi

./bin/dataset_QA -N -1 -picofile $input_list -outfile $out_name -system $system
	

