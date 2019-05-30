# 6.0002 Problem Set 5
# Graph optimization
# Name:
# Collaborators:
# Time:

import unittest

#
# A set of data structures to represent graphs
#

class Node(object):
    """Represents a node in the graph"""
    def __init__(self, name):
        self.name = str(name)

    def get_name(self):
        return self.name

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        # This function is necessary so that Nodes can be used as
        # keys in a dictionary, even though Nodes are mutable
        return self.name.__hash__()


class Edge(object):
    """Represents an edge in the dictionary. Includes a source and
    a destination."""
    def __init__(self, src, dest):
        self.src = src
        self.dest = dest

    def get_source(self):
        return self.src

    def get_destination(self):
        return self.dest

    def __str__(self):
        return '{}->{}'.format(self.src, self.dest)


class WeightedEdge(Edge):
    def __init__(self, src, dest, total_distance, outdoor_distance):
        self.src = src
        self.dest = dest
        self.total_distance = total_distance
        self.outdoor_distance = outdoor_distance

    def get_total_distance(self):
        return self.total_distance

    def get_outdoor_distance(self):
        return self.outdoor_distance

    def __str__(self):
        string = str(self.src) + "->" + str(self.dest) + " (" + str(self.total_distance) + ", " + str(self.outdoor_distance) + ")"
        return string


class Digraph(object):
    """Represents a directed graph of Node and Edge objects"""
    def __init__(self):
        self.nodes = set([])
        self.edges = {}  # must be a dict of Node -> list of edges
        self.neighbors = {} # dict of Node -> Neighbor Nodes

    def __str__(self):
        edge_strs = []
        for edges in self.edges.values():
            for edge in edges:
                edge_strs.append(str(edge))
        edge_strs = sorted(edge_strs)  # sort alphabetically
        return '\n'.join(edge_strs)  # concat edge_strs with "\n"s between them

    def get_edges_for_node(self, node):
        return self.edges[node]
    
    def childrenOf(self, node):
        return self.neighbors[node]

    def has_node(self, node):
        return node in self.nodes
    
    def isNodeValid(self, node):
        return node in self.edges

    def add_node(self, node):
        """Adds a Node object to the Digraph. Raises a ValueError if it is
        already in the graph."""
        if node.get_name() in self.edges:
            raise ValueError('Duplicate node')
        else:
            self.edges[node.get_name()] = []
            self.neighbors[node.get_name()] = []
            self.nodes.add(node)

    def add_edge(self, edge):
        """Adds an Edge or WeightedEdge instance to the Digraph. Raises a
        ValueError if either of the nodes associated with the edge is not
        in the  graph."""
        if not (type(edge) == Edge or type(edge) == WeightedEdge):
            raise ValueError('Not an edge!')
        src = edge.get_source()
        dest = edge.get_destination()
        if not (src.get_name() in self.edges and dest.get_name() in self.edges):
            raise ValueError('Node not in graph')
        self.neighbors[src.get_name()].append(dest.get_name()) #creates nieghborhood map
        self.edges[src.get_name()].append(edge)  #how the assignment wants it done

class Graph(Digraph):
    def add_edge(self, edge):
        """Adds an Edge or WeightedEdge instance to the Graph. Overrides
        the existing add_edge function in the Digraph object."""
        if not (type(edge) == Edge or type(edge) == WeightedEdge):
            raise ValueError('Not an edge!')
        Digraph.add_edge(self, edge)
        if type(edge) == Edge:
            rev = Edge(edge.get_destination(), edge.get_source())
        elif type(edge) == WeightedEdge:
            rev = WeightedEdge(edge.get_destination(), edge.get_source(), edge.get_total_distance(), edge.get_outdoor_distance())
        else:
            raise ValueError('Not an edge!')
        Digraph.add_edge(self, rev)

# ================================================================
# Begin tests -- you do not need to modify anything below this line
# ================================================================

class TestGraph(unittest.TestCase):

    def setUp(self):
        self.g = Digraph()
        self.na = Node('a')
        self.nb = Node('b')
        self.nc = Node('c')
        self.g.add_node(self.na)
        self.g.add_node(self.nb)
        self.g.add_node(self.nc)
        self.e1 = WeightedEdge(self.na, self.nb, 15, 10)
        self.e2 = WeightedEdge(self.na, self.nc, 14, 6)
        self.e3 = WeightedEdge(self.nb, self.nc, 3, 1)
        self.g.add_edge(self.e1)
        self.g.add_edge(self.e2)
        self.g.add_edge(self.e3)
        self.g2 = Graph()
        self.g2.add_node(self.na)
        self.g2.add_node(self.nb)
        self.g2.add_node(self.nc)
        self.g2.add_edge(self.e1)
        self.g2.add_edge(self.e2)
        self.g2.add_edge(self.e3)

    def test_weighted_edge_str(self):
        self.assertEqual(str(self.e1), "a->b (15, 10)")
        self.assertEqual(str(self.e2), "a->c (14, 6)")
        self.assertEqual(str(self.e3), "b->c (3, 1)")

    def test_weighted_edge_total_distance(self):
        self.assertEqual(self.e1.get_total_distance(), 15)
        self.assertEqual(self.e2.get_total_distance(), 14)
        self.assertEqual(self.e3.get_total_distance(), 3)

    def test_weighted_edge_outdoor_distance(self):
        self.assertEqual(self.e1.get_outdoor_distance(), 10)
        self.assertEqual(self.e2.get_outdoor_distance(), 6)
        self.assertEqual(self.e3.get_outdoor_distance(), 1)

    def test_add_edge_to_nonexistent_node_raises(self):
        node_not_in_graph = Node('q')
        no_src = WeightedEdge(self.nb, node_not_in_graph, 5, 5)
        no_dest = WeightedEdge(node_not_in_graph, self.na, 5, 5)

        with self.assertRaises(ValueError):
            self.g.add_edge(no_src)
        with self.assertRaises(ValueError):
            self.g.add_edge(no_dest)

    def test_add_existing_node_raises(self):
        with self.assertRaises(ValueError):
            self.g.add_node(self.na)

    def test_digraph_str(self):
        expected = "a->b (15, 10)\na->c (14, 6)\nb->c (3, 1)"
        self.assertEqual(str(self.g), expected)

    def test_graph_str(self):
        expected = "a->b (15, 10)\na->c (14, 6)\nb->a (15, 10)\nb->c (3, 1)\nc->a (14, 6)\nc->b (3, 1)"
        self.assertEqual(str(self.g2), expected)

    def test_add_non_edge_to_graph(self):
        with self.assertRaises(ValueError):
            self.g2.add_edge("Not an Edge")
        with self.assertRaises(ValueError):
            self.g.add_edge("ALSO Not an Edge")


if __name__ == "__main__":
    unittest.main()
