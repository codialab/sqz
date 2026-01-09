# D
- pop D{uv}
- remove D{uv} (i,a_0|a_1) (with left/right neighbor)
- add entry
- add entries
- get left neighbor
- get right neighbor

# Occurrence
- get self loops of a section

Pseudocode:

```python
for uv in D from highest freq to freq of 2:
  new meta-node q
  D_q <- pop D[uv]
  for each (i, a_1|a_2) in D_q:
    if left_neighbor[(u, i, a_1)] not v:
      (w, a_0) <- left_neighbor[(u, i, a_1)]
      remove D[wu] with (i, a_0|a_1)
      add D[wq] with (i, a_0|a_1)
  if u == v:
    Dvu <- self_loop_extract(D_q)
  else:
    Dvu <- pop D[vu]
  if vu in D or u == v:
    separate Dvu into Dqq, Dvq, Dqu, Dvu
    add Dqq as D[wq]
    add Dvq as D[vq]
    add Dqu as D[qu]
    add Dvu as D[vu]
  for each (i, a_1|a_2) in D_q:
    if right_neighbor[(v, i, a_2)] not u:
      (w, a_3) <- right_neighbor[(w, i, a_2)]
      remove D[vq] with (i, a_2|a_3)
      add D[wq] with (i, a_1|a_3)
  add (q, u, v, D_q) to Q
```
