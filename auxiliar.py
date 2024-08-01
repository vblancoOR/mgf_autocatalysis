# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 21:01:54 2024

@author: Trabajador
"""
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np


# =============================================================================
def convierteDiccionarioDeTiempo(xx):
    ts = list(set([j for i, j in xx]))
    ts.sort()
    dictXx = {}
    for t in ts:
        lista_t = []
        for i, j in xx:
            if j == t:
                lista_t.append(i)
        dictXx[t] = lista_t
    
    return dictXx
# =============================================================================


# =============================================================================
def convierteDiccionarioDeTiempoYNombre(xx, names):
    dictNames = dict(zip(range(len(names)), names))
    ts = list(set([j for i, j in xx]))
    ts.sort()
    dictXx = {}
    for t in ts:
        lista_t = []
        for i, j in xx:
            if j == t:
                lista_t.append(dictNames[i])
        dictXx[t] = lista_t
    
    return dictXx
# =============================================================================

# =============================================================================
def modifyFlow(dictFlow, names):
    dictNames = dict(zip(range(len(names)), names))
    keys = list(dictFlow.keys())
    values = list(dictFlow.values())
    keys = [(dictNames[i[0]], i[1]) for i in keys]
    sol = dict(zip(keys, values))
    return sol
# =============================================================================


# =============================================================================
def cambiarFormatoEspecies(ys):
    sol = [('s' + str(y[0] + 1), y[1]) for y in ys]
    return sol
# =============================================================================

# =============================================================================
def cambiarFormatoReacciones(zs, reacciones):

    correspondenciaDict = dict(zip(range(len(reacciones)), reacciones))

    sol = []
    for z in zs:
        reaccion = correspondenciaDict[z[0]]
        ayuda = (reaccion.split(' -> ')[0], reaccion.split(' -> ')[1], z[1])
        sol.append(ayuda)
    
    return sol
# =============================================================================


# =============================================================================
def creaGrafo(ys, zs):
    # Tiempos
    ts = range(max([i[1] for i in ys]) + 1)
    nombresSpecies = list(set([i[0] for i in ys]))
    nombresSpecies.sort()
    nombresReactives = list(set([i[0] for i in zs]))
    nombresReactives.sort()
    nombresProducts = list(set([i[1] for i in zs]))
    nombresProducts.sort()

    species = []
    reactivos = []
    productos = []
    for y in ys:
        species.append((y[0], y[1], 's'))
    for z in zs:
        reactivos.append((z[0], z[2], 'r'))
        productos.append((z[1], z[2], 'p'))
    species = list(set(species))
    reactivos = list(set(reactivos))
    productos = list(set(productos))

    G = nx.DiGraph()
    ps_anterior = []
    for t in ts:

        ns = [i for i in species 
              if i[1] == t]
        rs = [i for i in reactivos 
              if i[1] == t]
        ps = [i for i in productos 
              if i[1] == t]
        reacciones = [((z[0], z[2], 'r'),(z[1], z[2], 'p')) 
                      for z in zs 
                      if z[2] == t]
        for n in ns:
            G.add_node(n, typ = "specie")
        for r in rs:
            G.add_node(r, typ = "reactive")
        for p in ps:
            G.add_node(p, typ = "product")
        for i in reacciones:
            G.add_edge(i[0], i[1], typ = "reaction")
            
        for p in ps_anterior:
            for n in ns:
                if n[0] in p[0]:
                    G.add_edge(p, n, typ = "normal")
            
        ps_anterior = ps
        
        for n in ns:
            for r in rs:
                if n[0] in r[0]:
                    G.add_edge(n, r, typ = "normal")
            ps_anterior.append(n)

    return G
# =============================================================================



# =============================================================================
def reduceGraph(G, ys, zs):
    
    nombresSpecies = list(set([i[0] for i in ys]))
    nombresSpecies.sort()
    nombresReactives = list(set([i[0] for i in zs]))
    nombresReactives.sort()
    nombresProducts = list(set([i[1] for i in zs]))
    nombresProducts.sort()
    
    for s in nombresSpecies:
        listaAContraer = []
        for i in G.nodes:
            if G.nodes[i]["typ"] == "specie" and i[0] == s:
                listaAContraer.append(i)
        if len(listaAContraer) > 1:
            head, *tail = listaAContraer
            for j in tail:
                G = nx.contracted_nodes(G, head, j, self_loops=False)
            del G.nodes[head]['contraction']
    
    for r in nombresReactives:
        listaAContraer = []
        for i in G.nodes:
            # if i[0] == r and G.nodes[i]["typ"] == "reactive" or G.nodes[i]["typ"] == "product":
            if G.nodes[i]["typ"] == "reactive" and i[0] == r:
                listaAContraer.append(i)
        if len(listaAContraer) > 1:
            head, *tail = listaAContraer
            for j in tail:
                G = nx.contracted_nodes(G, head, j, self_loops=False)
            del G.nodes[head]['contraction']
            

    for p in nombresProducts:
        listaAContraer = []
        for i in G.nodes:
            # if i[0] == p and G.nodes[i]["typ"] == "reactive" or G.nodes[i]["typ"] == "product":
            if G.nodes[i]["typ"] == "product" and i[0] == p:
                listaAContraer.append(i)
        if len(listaAContraer) > 1:
            # print(listaAContraer)
            head, *tail = listaAContraer
            for j in tail:
                G = nx.contracted_nodes(G, head, j, self_loops=False)
            del G.nodes[head]['contraction']
            

    for i in G.edges:
        if 'contraction' in G.edges[i]:
            del G.edges[i]['contraction']
    
    return G
# =============================================================================


# =============================================================================
def filterSubGraphsByTypes(G, node_type1, node_type2):
    
    nodes = []
    for i in G.nodes:
        if (G.nodes[i]["typ"] == node_type1 or 
            G.nodes[i]["typ"] == node_type2):
            nodes.append(i)
            
    H = G.subgraph(nodes)
    subgraphsList = [list(w.nodes) 
                 for w in [H.subgraph(c) 
                           for c in
                           nx.connected_components(H.to_undirected())]]
    
    subgraphs = []
    for i in subgraphsList:
        subgraph = G.subgraph(i)
        subgraphs.append(subgraph)
    return subgraphs
# =============================================================================

# =============================================================================  
def modifyGraph(G, subgraphDict):
    for key, subgraph in subgraphDict.items():
        # Add contracted node representing the subgraph
        G.add_node(key, typ='contracted')
        for i in subgraph:
            # if i in G.nodes:
            G = nx.contracted_nodes(G, key, i, self_loops=False)
            # if 'contraction' in G.nodes[key]:
            del G.nodes[key]['contraction']

    for i in G.edges:
        if 'contraction' in G.edges[i]:
            del G.edges[i]['contraction']
            
    return G
# =============================================================================  
      

# =============================================================================  
def drawModify(G, subgraphDict):
    # Position for the nodes
    # -----------
    pos = graphviz_layout(G, prog="neato")
    # -----------

    # Labels for nodes
    # -----------
    labeldict = {}
    for i, d in G.nodes(data=True):
        if d['typ'] == "specie":
            labeldict[i] = i[0]
        elif d['typ'] == "contracted":
            labeldict[i] = ""
    # -----------

    # Size for nodes
    # -----------
    the_base_size = 350
    the_base_size_reaction = 10000
    sizeList = []
    for i, d in G.nodes(data=True):
        if d['typ'] == "specie":
            size = len(labeldict[i]) * the_base_size
        elif d['typ'] == "contracted":
            size = len(subgraphDict[i]) * the_base_size_reaction
        sizeList.append(size)
    # -----------

    # -----------
    # newPos = {}
    # scale_factor = 0.1  # Adjust this factor as needed
    # for idx, value in enumerate(pos.items()):
    #     size = sizeList[idx]
    #     newPos[value[0]] = (value[1][0] * size * scale_factor, value[1][1] * size * scale_factor)
    # pos = newPos
    # -----------

    nx.draw_networkx_nodes(G, pos, node_color= "none",
                           node_size = sizeList, 
                           node_shape = "o", 
                           alpha=1,
                           edgecolors = 'k',
                           linewidths=2)
    
    nx.draw_networkx_edges(G, pos, width=2.0, 
                           alpha=1, 
                           edge_color = "k", 
                           node_size = sizeList, 
                           node_shape = "o",
                           arrowsize= 20,
                           connectionstyle = "arc3")
    
    nx.draw_networkx_labels(G, pos, labeldict, alpha = 1,
                            font_weight = 700)

    
    for key, subgraph in subgraphDict.items():
        center = pos[key]
        radio = len(subgraphDict[i]) * the_base_size_reaction
        radio = radio/4000
        drawSubgraphInCircleV0(subgraph, center, radio)
# =============================================================================  


# =============================================================================  
def dibujaModificado(G, subgraphDict):
    # Position for the nodes
    # -----------
    pos = graphviz_layout(G, prog="neato")
    # -----------
   
    # Adjust node positions iteratively to prevent overlap
    # Draw edges
    # -----------
    # for edge in G.edges():
    #     start, end = adjustPositions(edge, pos, G, subgraphDict)
    #     arrow = FancyArrowPatch(start, end, arrowstyle='-|>', 
    #                             mutation_scale=40, color='k', zorder=3,
    #                             connectionstyle = "arc3", linewidth = 5,
    #                             alpha = 1)
    #     plt.gca().add_patch(arrow)

    for edge in G.edges():
        start, end = adjustPositions(edge, pos, G, subgraphDict)
        arrow = FancyArrowPatch(start, end, arrowstyle='-|>', 
                                mutation_scale=20, color='k', zorder=3,
                                connectionstyle = "arc3", linewidth = 2,
                                alpha = 1)
        plt.gca().add_patch(arrow)
    # -----------

    # Draw nodes
    # -----------
    # for i, d in G.nodes(data=True):
    #     if d["typ"] == 'specie':
    #         plt.text(pos[i][0], pos[i][1], str(i[0]), ha='center', va='center',
    #                   fontsize=30, fontweight='bold', zorder=3,
    #                   bbox=dict(boxstyle='round', facecolor='tomato',
    #                             edgecolor='black', alpha = 1, linewidth = 3))
    #         plt.scatter(*pos[i], s=0)
    #     elif d["typ"] == 'contracted':
    #         numberNodes = len(subgraphDict[i])
    #         radius = numberNodes * 10
    #         circle = plt.Circle(pos[i],
    #                             radius,
    #                             color = "royalblue",
    #                             linewidth=1,
    #                             alpha =0.5,
    #                             zorder=3, 
    #                             fill=True)
    #         plt.gca().add_patch(circle)

    for i, d in G.nodes(data=True):
        if d["typ"] == 'specie':
            plt.text(pos[i][0], pos[i][1], str(i[0]), ha='center', va='center',
                      fontweight='bold', zorder=3,
                      bbox=dict(boxstyle='round', facecolor='tomato',
                                edgecolor='black', alpha = 1, linewidth = 2))
            plt.scatter(*pos[i], s=0)
        elif d["typ"] == 'contracted':
            numberNodes = len(subgraphDict[i])
            radius = numberNodes * 5
            circle = plt.Circle(pos[i],
                                radius,
                                color = "royalblue",
                                linewidth=1,
                                alpha =0.5,
                                zorder=3, 
                                fill=True)
            plt.gca().add_patch(circle)
    # -----------
    
    
    # -----------
    # for key, subgraph in subgraphDict.items():
    #     center = pos[key]
    #     radio = len(subgraphDict[i]) * 10
    #     drawSubgraphInCircleV1(subgraph, center, radio)
    for key, subgraph in subgraphDict.items():
        center = pos[key]
        radio = len(subgraphDict[i]) * 5
        drawSubgraphInCircleV1(subgraph, center, radio)
    # -----------
    # Set axis off
    plt.axis('off')
    plt.show()
# ============================================================================= 



# ============================================================================= 
def adjustPositions(edge, pos, G, subgraphDict):
    segment = (pos[edge[0]], pos[edge[1]])
    if G.nodes[edge[0]]["typ"] == 'specie':
        # length_start = 5
        length_start = 2
    # elif G.nodes[edge[0]]["typ"] == 'contracted':
    #     length_start = len(subgraphDict[edge[0]])  * 10
    elif G.nodes[edge[0]]["typ"] == 'contracted':
        length_start = len(subgraphDict[edge[0]]) * 5  
    if G.nodes[edge[1]]["typ"] == 'specie':
        # length_end = 5
        length_end = 2
    # elif G.nodes[edge[1]]["typ"] == 'contracted':
    #     length_end = len(subgraphDict[edge[1]]) * 
    elif G.nodes[edge[1]]["typ"] == 'contracted':
        length_end = len(subgraphDict[edge[1]]) * 5
    start, end = shortenSegmentByScalar(segment, length_start, length_end)
    return start, end
# ============================================================================= 



# =============================================================================
def drawSubgraphInCircleV1(subgraph, center, radio):
    # Position for the nodes
    # -----------
    numberNodes = len(subgraph)
    pos = giveMeCoorSubgraph(subgraph, center, radio, 5)
    # -----------
    
    # Draw nodes
    # -----------
    # for i, d in subgraph.nodes(data=True):
    #     if d["typ"] == 'reactive':
    #         plt.text(pos[i][0], pos[i][1], str(i[0]), ha='center', va='center',
    #                   fontsize=30, fontweight='bold', zorder=3,
    #                   bbox=dict(boxstyle='round', facecolor='gold',
    #                                     edgecolor='k', alpha = 1, linewidth = 3))
    #         plt.scatter(*pos[i], s=0)
    #     elif d["typ"] == 'product':
    #         plt.text(pos[i][0], pos[i][1], str(i[0]), ha='center', va='center',
    #                   fontsize=30, fontweight='bold', zorder=3,
    #                   bbox=dict(boxstyle='round', facecolor='limegreen',
    #                                     edgecolor='k', alpha = 1, linewidth = 3))
    #     plt.scatter(*pos[i], s=0)
    
    for i, d in subgraph.nodes(data=True):
        if d["typ"] == 'reactive':
            plt.text(pos[i][0], pos[i][1], str(i[0]), ha='center', va='center',
                      fontweight='bold', zorder=3,
                      bbox=dict(boxstyle='round', facecolor='gold',
                                        edgecolor='k', alpha = 1, linewidth = 2))
            plt.scatter(*pos[i], s=0)
        elif d["typ"] == 'product':
            plt.text(pos[i][0], pos[i][1], str(i[0]), ha='center', va='center',
                      fontweight='bold', zorder=3,
                      bbox=dict(boxstyle='round', facecolor='limegreen',
                                        edgecolor='k', alpha = 1, linewidth = 2))
        plt.scatter(*pos[i], s=0)
    # -----------
    
    # Draw edges
    # -----------
    # for edge in subgraph.edges():
    #     segment = (pos[edge[0]], pos[edge[1]])
    #     # length_start = len(str(edge[0])) + 5
    #     # length_end = len(str(edge[1])) + 5
    #     length_start = len(str(edge[0]))/numberNodes - numberNodes + 5
    #     length_end = len(str(edge[1]))/numberNodes - numberNodes +4
    #     # print(length_start, length_end)
    #     start, end = shortenSegmentByScalar(segment, length_start, length_end)
    #     arrow = FancyArrowPatch(start, end, arrowstyle='-|>', 
    #                             mutation_scale=40, color='k',
    #                             zorder=3, linewidth = 5, alpha = 1, 
    #                             connectionstyle = "arc3")
    #     plt.gca().add_patch(arrow)
    
    for edge in subgraph.edges():
        segment = (pos[edge[0]], pos[edge[1]])
        # length_start = len(str(edge[0])) + 5
        # length_end = len(str(edge[1])) + 5
        length_start = len(str(edge[0]))/numberNodes - numberNodes
        length_end = len(str(edge[1]))/numberNodes - numberNodes
        # print(length_start, length_end)
        start, end = shortenSegmentByScalar(segment, length_start, length_end)
        arrow = FancyArrowPatch(start, end, arrowstyle='-|>', 
                                  color='k',
                                zorder=3, linewidth = 2, alpha = 1, 
                                connectionstyle = "arc3")
        plt.gca().add_patch(arrow)
    # -----------
# =============================================================================
  


# ============================================================================= 
def shortenSegmentByScalar(segment, length_start, length_end):
    """
    Shorten a segment from extremes by specified lengths.

    Parameters:
        segment: tuple of tuples, representing the segment ((x1, y1), (x2, y2))
        length_start: float, the length by which to shorten the segment
        length_end: float, the length by which to shorten the segment

    Returns:
        tuple of tuples, the shortened segment
    """
    (x1, y1), (x2, y2) = segment
    dx = x2 - x1
    dy = y2 - y1
    segment_length = (dx ** 2 + dy ** 2) ** 0.5

    # Calculate the new endpoint based on the chosen end to shorten from
    new_x1 = x1 + dx * (length_start / segment_length)
    new_y1 = y1 + dy * (length_start / segment_length)
    new_x2 = x2 - dx * (length_end / segment_length)
    new_y2 = y2 - dy * (length_end / segment_length)
    new_segment = ((new_x1, new_y1), (new_x2, new_y2))

    return new_segment
# ============================================================================= 


# =============================================================================
def drawSubgraphInCircleV0(subgraph, center, radio):
    
    map_color = []
    for i in subgraph:
        if i[2] == "r":
            color = "lightsteelblue"
        elif  i[2] == "p":
            color = "royalblue"
        map_color.append(color)
    
    pos = giveMeCoorSubgraph(subgraph, center, radio, 0.2)
    
    # Labels for nodes
    # -----------
    labeldict = {}
    for i in subgraph:
        labeldict[i] = i[0]
    # ----------- 
    
    # Size for nodes
    # -----------
    the_base_size = 350
    sizeList = []
    for i in subgraph:
        size = len(labeldict[i]) * the_base_size
        sizeList.append(size)
    # -----------
    
    GG = nx.DiGraph()
    for i in subgraph:
        GG.add_node(i)

    nx.draw_networkx_nodes(subgraph, pos, 
                           node_color= map_color,
                           node_size = sizeList,
                           node_shape = "o", 
                           alpha=0.7,
                           edgecolors = 'blue')
    
    nx.draw_networkx_edges(subgraph, pos, 
                           width=2.0, 
                           alpha=1, 
                           edge_color = "blue", 
                           node_size = sizeList, 
                           node_shape = "o",
                           arrowsize= 20,
                           connectionstyle = "arc3")
    
    nx.draw_networkx_labels(subgraph, pos, labeldict,
                            alpha = 1,
                            font_weight = 700)   
# =============================================================================


# =============================================================================
# Graph: graph
# center: center of the circle
# Radius: Radius of the circle
# Padding: Padding to keep nodes away from the circle border
def giveMeCoorSubgraph(SubgraphList, center, radius, padding):
    
    # Define parameters
    num_nodes = len(SubgraphList)

    # Calculate node positions
    # --
    # Evenly distribute nodes around the circle
    theta = np.linspace(0, 2*np.pi, num_nodes, endpoint=False)
    # Use a slightly smaller radius to keep nodes inside the circle
    r = radius - padding  
    positions = np.array([(center[0] + r * np.cos(t), center[1] + r * np.sin(t)) 
                          for t in theta])
    # --

    # Map node indices to positions
    node_positions = {node: positions[i] 
                      for i, node in enumerate(SubgraphList)}
    
    return node_positions
# =============================================================================




