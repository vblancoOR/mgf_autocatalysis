# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 21:18:31 2024

@author: Trabajador
"""

from networkx.drawing.nx_pydot import graphviz_layout
import networkx as nx
import auxiliar as aux


# =============================================================================
def dibujaAuxiliar(G):
    # Position for the nodes
    pos = graphviz_layout(G, prog="neato")
    # pos = nx.spring_layout(G)
    
    
    labeldict = {}
    for i in G.nodes:
        labeldict[i] = i[0]

    species = []
    reacciones = []
    for i in G.nodes():
        if G.nodes()[i]['typ'] == "specie":
            species.append(i)
        else:
            reacciones.append(i)
            
    the_base_size = 350
    sizeList = [len(labeldict[v]) * the_base_size for v in G.nodes()]
    
    # -----------
    # newPos = {}
    # scale_factor = 0.01  # Adjust this factor as needed
    # for idx, value in enumerate(pos.items()):
    #     size = sizeList[idx]
    #     newPos[value[0]] = (value[1][0] * size * scale_factor, value[1][1] * size * scale_factor)
    # pos = newPos
    # -----------

    colorMapNodes = []
    for i in G.nodes():
        if G.nodes()[i]['typ'] == "specie":
            colorMapNodes.append('red')
        elif G.nodes()[i]['typ'] == "reactive": 
            colorMapNodes.append('cornflowerblue')
        else:
            colorMapNodes.append('blue')

            
    colorMapEdges = []
    for i in G.edges():
        if G.edges()[i]['typ'] == "reaction":
            colorMapEdges.append('royalblue')
        else: 
            colorMapEdges.append('red')
            
            
    nx.draw_networkx_nodes(G, pos, node_color= colorMapNodes,
                            node_size = sizeList, node_shape = "o", alpha=0.2)
    
    nx.draw_networkx_edges(G, pos, width=2.0, 
                           alpha=0.6, 
                           edge_color = colorMapEdges, 
                           node_size = sizeList, 
                           node_shape = "o",
                           arrowsize= 15,
                           connectionstyle = "angle3")


    nx.draw_networkx_labels(G, pos, labeldict, alpha = 0.8,
                            font_weight = 700)
# =============================================================================





# =============================================================================
def dibujaAuxiliarPasos(G):
    # Position for the nodes
    pos = graphviz_layout(G, prog="dot")
    # pos = nx.spring_layout(G)
    
    labeldict = {}
    for i in G.nodes:
        labeldict[i] = i[0] + "\n" + "t: " +  str(i[1])
        
    ayudadict = {}
    for i in G.nodes:
        ayudadict[i] = i[0]
        

    species = []
    reacciones = []
    for i in G.nodes():
        if G.nodes()[i]['typ'] == "specie":
            species.append(i)
        else:
            reacciones.append(i)
            
    the_base_size = 350
    sizeDict = [len(ayudadict[v]) * the_base_size for v in G.nodes()]

    colorMapNodes = []
    for i in G.nodes():
        if G.nodes()[i]['typ'] == "specie":
            colorMapNodes.append('red')
        elif G.nodes()[i]['typ'] == "reactive": 
            colorMapNodes.append('cornflowerblue')
        else:
            colorMapNodes.append('blue')   
            
    colorMapEdges = []
    for i in G.edges():
        if G.edges()[i]['typ'] == "reaction":
            colorMapEdges.append('royalblue')
        else: 
            colorMapEdges.append('red')
            
            
    nx.draw_networkx_nodes(G, pos, node_color= colorMapNodes,
                            node_size = sizeDict, node_shape = 'o', alpha=0.2)

    nx.draw_networkx_edges(G, pos, width=2.0, 
                           alpha=0.6, 
                           edge_color = colorMapEdges, 
                           node_size = sizeDict, 
                           node_shape = "o",
                           arrowsize= 15,
                           connectionstyle = "angle3")

    nx.draw_networkx_labels(G, pos, labeldict, alpha = 0.8,
                            font_weight = 700 )
# =============================================================================











# =============================================================================  
def dibujaPasos(yValues, zValues):
    # Create Graph
    G = aux.creaGrafo(yValues, zValues)
    # Dibuja el grafo
    dibujaAuxiliarPasos(G)
# =============================================================================




# =============================================================================  
def dibuja(yValues, zValues):
    # Create Graph
    G = aux.creaGrafo(yValues, zValues)
    # Reduce graph
    G = aux.reduceGraph(G, yValues, zValues)
    # Dibuja el grafo
    dibujaAuxiliar(G)
# =============================================================================



# =============================================================================  
def dibujaRedDePetriV0(yValues, zValues):
    # Create Graph
    G = aux.creaGrafo(yValues, zValues)
    # Reduce graph
    G = aux.reduceGraph(G, yValues, zValues)
    # Got reactions subgraph
    subgraphs = aux.filterSubGraphsByTypes(G, "reactive", "product")
    # Dictionry sugraph
    subgraphDict = {}
    for i in range(len(subgraphs)):
        subgraphDict["A" + str(i)] = subgraphs[i]
    # Turn reactions into nodes
    G = aux.modifyGraph(G, subgraphDict)
    # Draw graph and subgrahs
    aux.drawModify(G, subgraphDict)
# =============================================================================  


# =============================================================================  
def dibujaRedDePetriV1(yValues, zValues):
    # Create Graph
    G = aux.creaGrafo(yValues, zValues)
    # Reduce graph
    G = aux.reduceGraph(G, yValues, zValues)
    # Got reactions subgraph
    subgraphs = aux.filterSubGraphsByTypes(G, "reactive", "product")
    # Dictionry sugraph
    subgraphDict = {}
    for i in range(len(subgraphs)):
        subgraphDict["A" + str(i)] = subgraphs[i]
    # Turn reactions into nodes
    G = aux.modifyGraph(G, subgraphDict)
    # Draw graph and subgrahs
    aux.dibujaModificado(G, subgraphDict)
# =============================================================================  
























