A Case for Work-stealing on FPGAs with OpenCL atomics 
==================================================

This repository contains reproducible code for the FPGA16 paper "A Case for Work-stealing on FPGAs with OpenCL Atomics" by Nadesh Ramanathan, Felix Winterstein, John Wickerson and George Constantinides. The key message of this work is to show that atomic operations can be used in a high-level synthesis(HLS) context to generate efficient custom hardware for irregular computations, K-means clustering in this case. Irregular workloads can be dynamic distributed via work-stealing and we implement work-stealing via OpenCL for FPGAs.

Code Overview   
-------------
Each OpenCL design (baseline /stealing) contains host and kernel code in "host/src" and "device" respectively. The host code takes input from the "golden_ref", where we have a number of datasets. For our paper, we only used one dataset which is 'n=1024*1024', 'k=128' and 'stv_dev=0.1' (evident in "main.cpp"). First, the host code generates either a balanced or unbalanced tree that is instantiated as OpenCL buffers and passed as OpenCL kernel argument. Then, the kernel runs 16 clustering iterations before writing its results into some OpenCL buffers that area subsequently read by the host code. Finally, the host code prints the final center results, execution time and the load distribution among OpenCL work-items. 

Interestingly, some form of load balancing is already obvious in simulation, if you just plan to have a high-level idea of this work and do not want to synthesize the entire kernel to understanding our work. Simply run ./build_sw.sh and compare the resulting outputs for the baseline and stealing designs. It must be warned that neither the simulation outputs nor the execution time potrait the actual load balancing effects, since we have induced "false/unwanted" barriers to trick the tool into doing some stealing during simulation (look for "SYN" preprocessing definition to have an idea about where the induced simulation barriers are). However, the induced barrier does give us an idea about what work-stealing can do to address workload imbalance. The actual load balancing results can be seen in our paper or after post-synthesis. Run ./build_hw.sh and the entire design will be compiled, synthesised and executed for you. 

Directory structure
-------------------
1) baseline: Contains the bline aseline implementation as described in the paper
2) work_stealing: Contains the work-stealing implementation as descirbed in the paper. This design uses atomics to load balancing via work-stealing. 
4) common: Contains Altera-specific OpenCL header files to allow compilation of OpenCL host and kernel code 

Experimentation Suggestions
-----------------------------
The parallelism degree in the current setup is 4 work-items. You can modify the definition of 'P' in host/src/filtering_algorithm_top.h for both the baseline and work-stealing implementations to vary the parallelism degree of your design. The node visited statistics will automatically output the new workload distribution (work done per work-item) for a different parallelism degree, even in simulation for the reasons described previously. We can fit up to 64 work-items for the baseline design and 32 work-items for the stealing design. You can run "make hw-estimate" to look at the estimated hardware resources (without going through synthesis) as you increase P.  

Contact
------------------
Please contact us at n_dot_ramanathan14_at_imperial_at_ic_ac_uk for any problems/bugs. 


(c) Nadesh Ramanathan, Imperial College London. The source code is distributed under a 3-clause BSD license (see LICENSE.txt). Please cite N. Ramanathan, F. Winterstein and G. Constantinides, "A Case for Work-stealing on FPGAs with OpenCL Atomics" in Proceedings of the 2016 ACM/SIGDA International Symposium on Field-Programmable Gate Arrays. 

