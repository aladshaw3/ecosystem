# 6.0002 Problem Set 5
# Graph optimization
# Name:
# Collaborators:
# Time:

#
# Finding shortest paths through MIT buildings
#
import unittest
from graph import Digraph, Node, WeightedEdge, Graph

#
# Problem 2: Building up the Campus Map
#
# Problem 2a: Designing your graph
#
# What do the graph's nodes represent in this problem? What
# do the graph's edges represent? Where are the distances
# represented?
#
# Answer:
#


# Problem 2b: Implementing load_map
def load_map(map_filename):
    """
    Parses the map file and constructs a directed graph

    Parameters:
        map_filename : name of the map file

    Assumes:
        Each entry in the map file consists of the following four positive
        integers, separated by a blank space:
            From To TotalDistance DistanceOutdoors
        e.g.
            32 76 54 23
        This entry would become an edge from 32 to 76.

    Returns:
        a Digraph representing the map
    """
    print("Loading map from file...")
    DG = Digraph()
    file = open(map_filename,"r")
    for line in file:
        line_list = line.split()
        if len(line_list) != 4:
            raise ValueError('Improper input file format')
        src = line_list[0]
        dest = line_list[1]
        tot_dist = 0.0
        out_dist = 0.0
        try:
            tot_dist = float(line_list[2])
        except:
            raise ValueError('Total Distance is not a number!')
        try:
            out_dist = float(line_list[3])
        except:
            raise ValueError('Outdoor Distance is not a number!')
        s_node = Node(src)
        d_node = Node(dest)
        cur_edge = WeightedEdge(s_node, d_node, tot_dist, out_dist)
        if not DG.has_node(s_node):
            DG.add_node(s_node)
        if not DG.has_node(d_node):
            DG.add_node(d_node)
        DG.add_edge(cur_edge)
    file.close()
    return DG

# Problem 2c: Testing load_map
# Include the lines used to test load_map below, but comment them out


#
# Problem 3: Finding the Shorest Path using Optimized Search Method
#
# Problem 3a: Objective function
#
# What is the objective function for this problem? What are the constraints?
#
# Answer:
#           Obj Func:   Minimize   { SUM( all edges, edge[i].tot_dist }
#
#           Contraint:  SUM( all edges, edge[i].out_dist) <= Max_Outdoor_Dist

def printPath(path):
    """Assumes path is a list of nodes"""
    result = ''
    for i in range(len(path)):
        result = result + str(path[i])
        if i != len(path) - 1:
            result = result + '->'
    return result

#Can be modified to consider the weights of each path
def DFS(graph, start, end, path, shortest, toPrint = False):
    """Assumes graph is a Digraph; start and end are nodes;
        path and shortest are lists of nodes
        Returns a shortest path from start to end in graph"""
    path = path + [start]
    if toPrint:
        print('Current DFS path:', printPath(path))
    if start == end:
        return path
    for node in graph.childrenOf(start):
        if node not in path: #avoid cycles
            if shortest == None or len(path) < len(shortest):
                newPath = DFS(graph, node, end, path, shortest,
                              toPrint)
                if newPath != None:
                    if toPrint:
                        print('newPath:', printPath(newPath))
                    shortest = newPath
        elif toPrint:
            print('Already visited', node)
    return shortest

def shortestPath(graph, start, end, toPrint = False):
    """Assumes graph is a Digraph; start and end are nodes
        Returns a shortest path from start to end in graph"""
    return DFS(graph, start, end, [], None, toPrint)


# Problem 3b: Implement get_best_path
def get_best_path(digraph, start, end, path, max_total_dist, max_dist_outdoors, best_dist, best_out, best_path):
    """
    Finds the shortest path between buildings subject to constraints.

    Parameters:
        digraph: Digraph instance
            The graph on which to carry out the search
        start: string
            Building number at which to start
        end: string
            Building number at which to end
        path: list composed of [[list of strings], int, int]
            Represents the current path of nodes being traversed. Contains
            a list of node names, total distance traveled, and total
            distance outdoors.
        max_total_dist: int
            Maximum total distance on a path
        max_dist_outdoors: int
            Maximum distance spent outdoors on a path
        best_dist: int
            The smallest distance between the original start and end node
            for the initial problem that you are trying to solve
        best_out: int
            The smallest distance between the original start and end node
            for outdoor distance traveled on the path
        best_path: list of strings
            The shortest path found so far between the original start
            and end node.

    Returns:
        A tuple with the shortest-path from start to end, represented by
        a list of building numbers (in strings), [n_1, n_2, ..., n_k],
        where there exists an edge from n_i to n_(i+1) in digraph,
        for all 1 <= i < k and the distance of that path.

        If there exists no path that satisfies max_total_dist and
        max_dist_outdoors constraints, then return None.
    """
    if digraph.isNodeValid(start) == True:
        path[0] += [start]
    else:
        return (None, 0, 0)
    # Compute length of path
    if len(path[0]) > 1:
        #dist between start and path[0][-2]
        for edge in digraph.get_edges_for_node(path[0][-2]):
            if edge.get_destination().get_name() == start:
                path[1] += edge.get_total_distance()
                path[2] += edge.get_outdoor_distance()

    if path[1] > max_total_dist or path[2] > max_dist_outdoors:
        return (None, 0, 0)

    result = (None, 0, 0)

    #Check for finished path
    if start == end:
        result = path
    else:
        #Iterate through the edges of the current node
        for edge in digraph.get_edges_for_node(start):
            #If current destination is unique, then go down that path
            nextNode = edge.get_destination().get_name()
            if nextNode not in path[0]:
                newPath = (None, 0, 0)
                #NOTE: Only make a search if path is strictly shorter than the best path so far
                if best_path == None or path[1] < best_dist:
                    newPath = get_best_path(digraph, nextNode, end, path, max_total_dist, max_dist_outdoors, best_dist, best_out, best_path)
                    
                    if newPath[0] == None:
                        if len(path[0]) > 1:
                            #dist between start and path[0][-2]
                            for e in digraph.get_edges_for_node(path[0][-2]):
                                if e.get_destination().get_name() == nextNode:
                                    path[1] = path[1] - e.get_total_distance()
                                    path[2] = path[2] - e.get_outdoor_distance()
                        if len(path[0]) > 0:
                            path[0].pop()
                    else:
                        if best_path == None:
                            best_path = newPath[0].copy()
                            best_dist = newPath[1]
                            best_out = newPath[2]
                            result = (best_path, best_dist, best_out)
                        else:
                            #NOTE: Only update if path is strictly shorter than the best path so far
                            if newPath[1] < best_dist:
                                best_path = newPath[0].copy()
                                best_dist = newPath[1]
                                best_out = newPath[2]
                            result = (best_path, best_dist, best_out)

                #Now have a valid result, but not necessarily the optimal result
                if best_path != None and newPath[0] != None:
                    if len(path[0]) > 1:
                        #dist between start and path[0][-2]
                        for e in digraph.get_edges_for_node(path[0][-2]):
                            if e.get_destination().get_name() == nextNode:
                                path[1] = path[1] - e.get_total_distance()
                                path[2] = path[2] - e.get_outdoor_distance()
                    if len(path[0]) > 0:
                        path[0].pop()
        #End edge loop
    return result


# Problem 3c: Implement directed_dfs
def directed_dfs(digraph, start, end, max_total_dist, max_dist_outdoors):
    """
    Finds the shortest path from start to end using a directed depth-first
    search. The total distance traveled on the path must not
    exceed max_total_dist, and the distance spent outdoors on this path must
    not exceed max_dist_outdoors.

    Parameters:
        digraph: Digraph instance
            The graph on which to carry out the search
        start: string
            Building number at which to start
        end: string
            Building number at which to end
        max_total_dist: int
            Maximum total distance on a path
        max_dist_outdoors: int
            Maximum distance spent outdoors on a path

    Returns:
        The shortest-path from start to end, represented by
        a list of building numbers (in strings), [n_1, n_2, ..., n_k],
        where there exists an edge from n_i to n_(i+1) in digraph,
        for all 1 <= i < k

        If there exists no path that satisfies max_total_dist and
        max_dist_outdoors constraints, then raises a ValueError.
    """
    result = get_best_path(digraph, start, end, [[],0,0], max_total_dist, max_dist_outdoors, 0, 0, None)
    if result[0] == None:
        raise ValueError('No result found!')
    if result[1] > max_total_dist:
        raise ValueError('Exceeds Maximum Total Distance!')
    if result[2] > max_dist_outdoors:
        raise ValueError('Exceeds Maximum Outdoor Distance!')
    return result[0]
    #return result

# My tests

#graph = load_map("mit_map.txt")
#graph = load_map("my_map.txt")
#print(graph)
#print()
#result = get_best_path(graph, 'a', 'g', [[],0,0], 100, 100, 0, 0, None) #Good
#result = directed_dfs(graph, 'a', 'g', 300, 0)

#print(graph)
#result = directed_dfs(graph, '10', '32', 100, 1000)
#result = get_best_path(graph, '3', '4', [['2'],0,0], 100, 100, 0, 0, None)
#result = shortestPath(graph, '2', '9', True )
#print("result")
#print(result)


# ================================================================
# Begin tests -- you do not need to modify anything below this line
# ================================================================

class Ps2Test(unittest.TestCase):
    LARGE_DIST = 99999

    def setUp(self):
        self.graph = load_map("mit_map.txt")

    def test_load_map_basic(self):
        self.assertTrue(isinstance(self.graph, Digraph))
        self.assertEqual(len(self.graph.nodes), 37)
        all_edges = []
        for _, edges in self.graph.edges.items():
            all_edges += edges  # edges must be dict of node -> list of edges
        all_edges = set(all_edges)
        self.assertEqual(len(all_edges), 129)


    def _print_path_description(self, start, end, total_dist, outdoor_dist):
        constraint = ""
        if outdoor_dist != Ps2Test.LARGE_DIST:
            constraint = "without walking more than {}m outdoors".format(
                outdoor_dist)
        if total_dist != Ps2Test.LARGE_DIST:
            if constraint:
                constraint += ' or {}m total'.format(total_dist)
            else:
                constraint = "without walking more than {}m total".format(
                    total_dist)

        print("------------------------")
        print("Shortest path from Building {} to {} {}".format(
            start, end, constraint))

    def _test_path(self,
                   expectedPath,
                   total_dist=LARGE_DIST,
                   outdoor_dist=LARGE_DIST):
        start, end = expectedPath[0], expectedPath[-1]
        self._print_path_description(start, end, total_dist, outdoor_dist)
        dfsPath = directed_dfs(self.graph, start, end, total_dist, outdoor_dist)
        print("Expected: ", expectedPath)
        print("DFS: ", dfsPath)
        self.assertEqual(expectedPath, dfsPath)

    def _test_impossible_path(self,
                              start,
                              end,
                              total_dist=LARGE_DIST,
                              outdoor_dist=LARGE_DIST):
        self._print_path_description(start, end, total_dist, outdoor_dist)
        with self.assertRaises(ValueError):
            directed_dfs(self.graph, start, end, total_dist, outdoor_dist)

    def test_path_one_step(self):
        self._test_path(expectedPath=['32', '56'])

    def test_path_no_outdoors(self):
        self._test_path(
            expectedPath=['32', '36', '26', '16', '56'], outdoor_dist=0)

    def test_path_multi_step(self):
        self._test_path(expectedPath=['2', '3', '7', '9'])

    def test_path_multi_step_no_outdoors(self):
        self._test_path(
            expectedPath=['2', '4', '10', '13', '9'], outdoor_dist=0)

    def test_path_multi_step2(self):
        self._test_path(expectedPath=['1', '4', '12', '32'])

    def test_path_multi_step_no_outdoors2(self):
        self._test_path(
            expectedPath=['1', '3', '10', '4', '12', '24', '34', '36', '32'],
            outdoor_dist=0)

    def test_impossible_path1(self):
        self._test_impossible_path('8', '50', outdoor_dist=0)

    def test_impossible_path2(self):
        self._test_impossible_path('10', '32', total_dist=100)

if __name__ == "__main__":
    unittest.main()
