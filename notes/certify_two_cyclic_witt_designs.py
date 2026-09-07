#!/usr/bin/env sage-python
"""Verify the two cyclic Witt designs and their degree-23 Mathieu groups."""
from sage.all import GF, PolynomialRing, Graph, PermutationGroup, SymmetricGroup
from sage.libs.gap.libgap import libgap

B=PolynomialRing(GF(2),"x");x=B.gen()
g=x**11+x**9+x**7+x**6+x**5+x+1
blocks=[]
for mask in range(1<<12):
    f=g*sum(x**j for j in range(12) if (mask>>j)&1)
    support=tuple(j for j in range(23) if f[j])
    if len(support)==7:blocks.append(support)
assert len(blocks)==253
graph=Graph()
graph.add_vertices(range(276))
graph.add_edges((point,23+number) for number,block in enumerate(blocks) for point in block)
automorphisms=graph.automorphism_group(partition=[list(range(23)),list(range(23,276))])
assert automorphisms.order()==10200960
group=PermutationGroup([[permutation(j)+1 for j in range(23)] for permutation in automorphisms.gens()])
assert group.order()==10200960 and group.is_simple()
assert group.is_transitive() and libgap.Transitivity(group)==4
ambient=SymmetricGroup(23)
translation=ambient(tuple(range(1,24)))
P=ambient.subgroup([translation])
assert translation in group
assert group.normalizer(group.subgroup([group(translation)])).order()==253
normalizer=ambient.normalizer(P)
assert normalizer.order()==506
assert ambient.normalizer(group)==ambient.subgroup([ambient(g) for g in group.gens()])
negation=ambient([(-j)%23+1 for j in range(23)])
other=ambient.subgroup([negation**(-1)*ambient(g)*negation for g in group.gens()])
assert other!=group and translation in other
assert group.intersection(other).order()==253
assert ambient.subgroup([*[ambient(g) for g in group.gens()],*other.gens()]).order()==ambient.order()//2
print("WITT_DESIGN_AUTOMORPHISM_GROUP: simple, degree23, order10200960, transitivity4")
print("NORMALIZERS: N_M23(P)=253, N_S23(P)=506, N_S23(M23)=M23")
print("TWO_COPIES: number containing a fixed P is 506/253=2")
print("NEGATION_EXCHANGES_COPIES: intersection23:11; generated groupA23")
print("PASS_TWO_CYCLIC_WITT_DESIGNS")
