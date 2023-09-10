using Metis
using SparseArrays 
using LinearAlgebra
using PyCall

rows=[1,2,3,4,6,5,5,7,8,9,10,11,12,13,14,15]
cols=[2,3,4,6,1,6,7,8,9,10,5,12,13,14,15,11]
vals=[1 for _=1:lastindex(rows)]
S=sparse(rows,cols,vals)
S=S+S'
G=Metis.graph(S)

ans=Metis.partition(G,3)
print(ans)



