ctl_dir=/home/io500/sys_structure_lab/

# 创建与perf通信的管道
ctl_fifo=${ctl_dir}perf_ctl.fifo
test -p ${ctl_fifo} && unlink ${ctl_fifo}
mkfifo ${ctl_fifo}
exec {ctl_fd}<>${ctl_fifo}
export CTL_FD=${ctl_fd} 

# 36 * 36
perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses     \
          --control fd:${ctl_fd}               \
          -o ./info1.txt ./bin/matrix_multiple 32 s

perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses     \
          --control fd:${ctl_fd}               \
          -o ./info1.txt --append ./bin/matrix_multiple 32 r

perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses     \
          --control fd:${ctl_fd}               \
          -o ./info1.txt --append ./bin/matrix_multiple 32 b 8

perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses     \
          --control fd:${ctl_fd}               \
          -o ./info1.txt --append ./bin/matrix_multiple 32 v   

perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses      \
          --control fd:${ctl_fd}               \
          -o ./info1.txt --append ./bin/matrix_multiple 32 w 8


# 4096 * 4096
perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses     \
          --control fd:${ctl_fd}               \
          -o ./info1.txt --append ./bin/matrix_multiple 4096 s

perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses     \
          --control fd:${ctl_fd}               \
          -o ./info1.txt --append ./bin/matrix_multiple 4096 r

perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses     \
          --control fd:${ctl_fd}               \
          -o ./info1.txt --append ./bin/matrix_multiple 4096 b 256

perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses     \
          --control fd:${ctl_fd}               \
          -o ./info1.txt --append ./bin/matrix_multiple 4096 v   

perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses     \
          --control fd:${ctl_fd}               \
          -o ./info1.txt --append ./bin/matrix_multiple 4096 v o

perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses      \
          --control fd:${ctl_fd}               \
          -o ./info1.txt --append ./bin/matrix_multiple 4096 w 256

perf stat -D -1 -e cpu-cycles,instructions,branches,branch-misses,cache-references,cache-misses      \
          --control fd:${ctl_fd}               \
          -o ./info1.txt --append ./bin/matrix_multiple 4096 w 256 o

# cleanup
perf_pid=$!
kill $perf_pid
rm -rf ${ctl_fifo}
unset CTL_FD

