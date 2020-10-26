# SymBi
Symmetric Continuous Subgraph Matching with Bidirectional Dynamic Programming

## Usages
./symbi <initial_data_graph_file> <graph_update_stream_file> <query_graph_file> [num_operations]

## Experiments in the paper
- Each script file in `scripts` directory reproduces the results of the experiments in the paper.
  - exp1_netflow_varying_query_size.sh: Figure 4
  - exp2_lsbench_varying_query_size.sh: Figure 5
  - exp3_netflow_varying_deletion_rate.sh: Figure 7
  - exp4_lsbench_varying_deletion_rate.sh: Figure 9
  - exp5_netflow_varying_insertion_rate.sh: Figure 10
  - exp6_lsbench_varying_insertion_rate.sh: Figure 12
  - exp7_lsbench_varying_dataset_size.sh: Figure 13
- The results will be stored in `results/<exp_name>` directory.

## Datasets
Put extracted zip file into `datasets` directory.
- Netflow ([download](https://drive.google.com/file/d/1g1NVK29k27V76JrSAVhKOZ2qvyP_kj--/view?usp=sharing))
- LSBench_x1 ([download](https://drive.google.com/file/d/1MmfMTE2SKPOJaaSibPC-5OPXh7brDulf/view?usp=sharing))
- LSBench_x5 ([download](https://drive.google.com/file/d/1k0eQdZYElAsBqk5NfhW414FsZuhDNsqG/view?usp=sharing))
- LSBench_x25 ([download](https://drive.google.com/file/d/1c_pGxT18H-BOL0cirdTLJ-0qbqcaDADx/view?usp=sharing))
