'''
Michael Norris
Take data from http://www.wormatlas.org/neuronalwiring.html and convert it to
networkx format.
'''

from xlrd import open_workbook
import networkx

#b1 = open_workbook('data/NeuronLineage_Part1.xls', on_demand=True)
#b2 = open_workbook('data/NeuronLineage_Part2.xls', on_demand=True)
b1 = open_workbook('data/NeuronConnect.xls', on_demand=True)
#b1 = open_workbook('data/NeuronFixedPoints.xls', on_demand=True)

graph = networkx.Graph()

#spreadsheets have 3 columns
#Neuron 1, Neuron 2, Relatedness

for book in [b1]:#[b1, b2]:
    for sheet in book.sheet_names():
        sheet = book.sheet_by_name(sheet)
        
        neuron1 = sheet.col(0)[1:]
        neuron2 = sheet.col(1)[1:]
        #relatedness = sheet.col(2)[2:]
        _type = sheet.col(2)[2:]
        nbr = sheet.col(3)[2:]

        for n1,n2,t,n in zip(neuron1, neuron2, _type, nbr):
            #graph.add_edge(n1.value, n2.value, relatedness=int(rel.value))
            graph.add_edge(n1.value, n2.value, 
                           type=t, nbr=n)

networkx.write_adjlist(graph, 'data/wormatlas.graphml')
