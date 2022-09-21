# MCQA: Muti-Constraint Qubit Allocation for Near-Term FTQC (developed in C++)

## Contact
This code is written by [Sunghye Park](shpark96@postech.ac.kr) and [Dohun Kim](dohunkim@postech.ac.kr)

If you have any questions feel free to contact us using shpark96@postech.ac.kr

(CSDL Lab. in POSTECH, South Korea. http://csdl.postech.ac.kr)

## Related paper
- Title: MCQA: Multi-Constraint Qubit Allocation for Near-FTQC Device
- Authors: S.Park, D.Kim, J.-Y.Sim, and S.Kang
- Conference: IEEE/ACM International Conference on Computer-Aided Design (ICCAD), 2022, to appear

## Overview
In this work, we present an effective heuristic approach for qubit allocation that considers not only the hardware topology but also other constraints for near-fault-tolerant quantum computing (near-FTQC).

## Usage

### Installation
```
[~/MQCA] mkdir build && cd build
[~/MCQA/build] cmake ..
[~/MCQA/build] make -j
```

### Run
```
[~/MCQA/tool] python3 run.py [circuit_index]
```
or
```
[~/MCQA/tool] python3 run.py all
```


If you want to know more about "run.py", put this command `python3 run.py help'

### System Requirements

* CMake >= 3.20
* Boost >= 1.60
* Sparsehash

## Files
| File      | Description |
| ----------- | ----------- |
| main.cpp      |  |
| circuit.h   |  |
| initial_mapping.cpp   |  |
| mapping.cpp   |  |
| mappingPreparation.cpp   |  |
| mappingFunction.cpp   |  |
| couplinggraph.cpp  |  |
| graph.cpp   |  |
| graphFunction.cpp   |  |
| layout.cpp   |  |
| latency.cpp   |  |
| outputwriter.cpp   |  |
| parser.cpp  |  |
| qasm-tools/QASM*   |  |
| ../examples/| Benchmark circuits for qubit allocation |


## Reference
If you use out mapping algorithm for your research, we would be thankful if you referred to it by citing the following publication
```
{
  title={MCQA: Multi-Constraint Qubit Allocation for Near-FTQC Device},
  author={Park, Sunghye and Kim, Dohun and Sim, Jae-Yoon and Kang, Seokhyeong},
  booktitle={Proceedings of the 41th IEEE/ACM International Conference on Computer-Aided Design},
  year={2022}
}
```
