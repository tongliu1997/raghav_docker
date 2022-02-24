#!/bin/bash
#SBATCH  --job-name=isobar_docker
#SBATCH  --array=0-19
#SBATCH  -o out-slurm/%a

#SBATCH   --partition=scavenge
#SBATCH   --time=01:00:00
#SBATCH   --mem=2G

n_events=-1

system_list=("Ru" "Zr")
system_id=$((  ${SLURM_ARRAY_TASK_ID}%2   ))
list_id=$(( ${SLURM_ARRAY_TASK_ID}/2  ))
system=${system_list[${system_id}]}

inp_list=./lists/${system}_list_${list_id}.list
while read line; do
    echo $line
done < ${inp_list}

output_dir=./out_root
output_name=${output_dir}/${system}_output_${list_id}.root
#singularity exec -B /gpfs/loomis/home.grace/tl543/Isobar_Data star_star.simg bash runjob.sh testfile.list test_output.root
#singularity exec -B /gpfs/loomis/project/caines/tl543/Isobar_Data star_star.simg bash runjob.sh ${inp_list} ${output_name} ${system}
