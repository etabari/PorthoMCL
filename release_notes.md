# PorthoMCL
Parallel implementation of OrthoMCL


## 1. OrthoMCL Bug

OrthoMCL's definition for pairs from [their documentation](https://docs.google.com/document/d/1RB-SqCjBmcpNq-YbOYdFxotHGuU7RK_wqxqDAMjyP_w/pub) is as follow:

```
- Input: BLAST results
- QiH(Ax,Ay) iff:
  - Ax != Ay
  - H(Ax,Ay) > Cutoff
  - For all taxa T != A and protein n in T:
     - H(Ax,Ay) >= H(Ax,Tn)
- IP(Ax,Ay) iff:
   - QiH(Ax,Ay) && QiH(Ay,Ax)
 ```