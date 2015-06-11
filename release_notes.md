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

where 

- `Ax`: protein `x` in taxa `A` 
- `Ay`: protein `y` in taxa `A` 
- `QiH(Ax,Ay)`: qualifying in-species hit (single way better hit) (br.csv in our implementation) 

Having said that the SQL command that supplies `QiH(Ax,Ay)` is as follows:

```sql
create table BetterHit  as
select s.query_id, s.subject_id,
       s.query_taxon_id as taxon_id,
       s.evalue_exp, s.evalue_mant
from SimilarSequences s, BestInterTaxonScore bis
where s.query_id != s.subject_id 
  and s.query_taxon_id = s.subject_taxon_id
  and s.query_id = bis.query_id
  and s.evalue_exp <= -5
  and s.percent_match >= 50
  and (s.evalue_mant < 0.001
       or s.evalue_exp < bis.evalue_exp
       or (s.evalue_exp = bis.evalue_exp and s.evalue_mant <= bis.evalue_mant))
-- . . . or Similarity for a protein with no BestInterTaxonScore
--       (i.e. an intrataxon match for a protein with no intertaxon
--        match in the database)
union
select s.query_id, s.subject_id, s.query_taxon_id as taxon_id, s.evalue_exp, s.evalue_mant
from SimilarSequences s
where s.query_taxon_id = s.subject_taxon_id 
  and s.evalue_exp <= -5
  and s.percent_match >= 50
  and s.query_id in 
     (SELECT distinct ust.query_id
      from UniqSimSeqsQueryId ust
      LEFT OUTER JOIN BestInterTaxonScore bis ON bis.query_id = ust.query_id
      WHERE bis.query_id IS NULL)
```