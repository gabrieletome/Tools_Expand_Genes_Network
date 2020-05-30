import numpy as np
import networkx
import rpack

#Return the dictionary with the layout of the graph
#Input: networkx.Graph, networkx type layout
def layout_many_components(graph, component_layout_func=networkx.layout.spring_layout, pad_x=1., pad_y=1.):
    components = _get_components_sorted_by_size(graph)
    component_sizes = [len(component) for component in components]
    bboxes = _get_component_bboxes(component_sizes, pad_x, pad_y)

    pos = dict()
    for component, bbox in zip(components, bboxes):
        component_pos = _layout_component(component, bbox, component_layout_func)
        pos.update(component_pos)

    return pos

#Sort component by nuomber of nodes
def _get_components_sorted_by_size(g):
    subgraphs = list(g.subgraph(c) for c in networkx.connected_components(g))
    return sorted(subgraphs, key=len)

def _get_component_bboxes(component_sizes, pad_x=1., pad_y=1.):
    dimensions = [_get_bbox_dimensions(n, power=0.8) for n in component_sizes]
    # rpack only works on integers; sizes should be in descending order
    dimensions = [(int(width + pad_x), int(height + pad_y)) for (width, height) in dimensions[::-1]]
    origins = rpack.pack(dimensions)
    bboxes = [(x, y, width-pad_x, height-pad_y) for (x,y), (width, height) in zip(origins, dimensions)]
    return bboxes[::-1]

def _get_bbox_dimensions(n, power=0.5):
    return (n**power, n**power)

#Return the position of nodes of a component
def _layout_component(component, bbox, component_layout_func):
    pos = component_layout_func(component)
    rescaled_pos = _rescale_layout(pos, bbox)
    return rescaled_pos

#Rescale layout based on the components
def _rescale_layout(pos, bbox):
    min_x, min_y = np.min([v for v in pos.values()], axis=0)
    max_x, max_y = np.max([v for v in pos.values()], axis=0)

    if not min_x == max_x:
        delta_x = max_x - min_x
    else: # graph probably only has a single node
        delta_x = 1.

    if not min_y == max_y:
        delta_y = max_y - min_y
    else: # graph probably only has a single node
        delta_y = 1.

    new_min_x, new_min_y, new_delta_x, new_delta_y = bbox
    new_pos = dict()
    for node, (x, y) in pos.items():
        new_x = (x - min_x) / delta_x * new_delta_x + new_min_x
        new_y = (y - min_y) / delta_y * new_delta_y + new_min_y
        new_pos[node] = (new_x, new_y)

    return new_pos
