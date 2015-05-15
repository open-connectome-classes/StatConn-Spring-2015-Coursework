from opendataset import graph
from blocks import blocks
for b in blocks:
    for n in blocks[b]:
        assert n in graph.edge.keys()
