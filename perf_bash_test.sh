ctl_dir=/home/io500/sys_structure_lab/

ctl_fifo=${ctl_dir}perf_ctl.fifo
test -p ${ctl_fifo} && unlink ${ctl_fifo}
mkfifo ${ctl_fifo}
exec {ctl_fd}<>${ctl_fifo}
echo ${ctl_fd}
export CTL_FD=${ctl_fd} 

perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses      \
          --control fd:${ctl_fd}               \
          -o ./info2.txt ./bin/matrix_multiple 4096 s

perf_pid=$!
kill $perf_pid
rm -rf ${ctl_fifo}
unset CTL_FD

