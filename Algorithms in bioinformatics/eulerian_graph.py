#!/usr/bin/env python3

"""
Author: Thierry Haddad
Student nr: 
Script to: 
"""
# Import statements
import random
from collections import deque

# Implement your functions here
def assert_dict(graph):
    if isinstance(graph, dict):
        return True
    else:
        print("Graph is not a dictionary type!")
        return False

def get_vertices(graph):
    return graph.keys()

def get_edges(graph):
    return graph.values()

def get_vertices_and_edges(graph):
    return graph.items()

def get_balance(graph):
    # {vertex: [#in, #out]}
    balance_dict = {}

    # Make new dict entry with known out's
    for vertex in get_vertices(graph):
        balance_dict[vertex] = [0, len(graph[vertex])]

    # Go through all edges to count in's
    for vertex, edges in get_vertices_and_edges(graph):
        for edge in edges:
            if edge in balance_dict:
                balance_dict[edge][0] += 1
            else:
                balance_dict[edge] = [1, 0]  # 1 in, 0 out
    return balance_dict

def has_eulerian_path(graph, balance_dict):
    # Every edge once
    all_edges = get_edges(balance_dict)
    edges = []
    # Find non-balanced vertices
    for edge in all_edges:
        if edge[0]!=edge[1]:
            edges.append(edge)
    # No un-balanced -> return True
    if len(edges) == 0:
        return True
    # Check if they are semi-balanced
    elif len(edges) <= 2:
        for edge in edges:
            # vertex(|ins-outs|) should equal 1
            if abs(edge[0]-edge[1]) != 1:
                return False
        # Return true if 0, 1 or 2 vertices are semi-balanced, rest balanced.
        return True
    # >2 un-balanced -> return False
    else:
        return False
    return False

def is_balanced(balance_dict):
    # Check balance
    for edges in get_edges(balance_dict):
        if not edges[0]==edges[1]:
            return False
    return True

def is_connected(graph):
    # check if all vertices are connected; no disconnect
    #TODO
    return True

def is_eulerian(balance_dict):
    if is_balanced(balance_dict):
        return True
    return False

def find_eulerian_cycle(graph, balance_dict):
    # Start/end same vertice
    # make subcycles and merge them
    if is_balanced(balance_dict):
        cycles = []
        copy_graph = graph.copy()
        while copy_graph != {}:
            path = []
            # Find random start vertex v
            start_vertex = random.choice(list(get_vertices(copy_graph)))
            # Find initial edge
            vertex = start_vertex
            path.append(vertex)
            choices = copy_graph[vertex]
            original_vertex = vertex
            vertex = random.choice(choices)
            choices.remove(vertex)
            copy_graph[original_vertex] = choices
            # Iterate untill current vertex == start vertex
            while vertex != start_vertex:
                path.append(vertex)
                original_vertex = vertex
                choices = copy_graph[vertex]
                vertex = random.choice(choices)
                choices.remove(vertex)
                if choices == []:
                    del copy_graph[original_vertex]
                else:
                    copy_graph[original_vertex] = choices
            cycles.append(path)
            vertices = list(get_vertices(copy_graph))
            for vertex in vertices:
                if copy_graph[vertex] == []:
                    del copy_graph[vertex]
        print(cycles)
        merge_paths(cycles)
        if len(cycles) >= 1:
            return True
        else:
            return False
    else:
        return False

def find_eulerian_path(graph, balance_dict):
    # Try to find way to visit every vertice, using every edge only once
    if not find_eulerian_cycle(graph, balance_dict):
        if has_eulerian_path(graph, balance_dict):
            paths = []
            copy_graph = graph.copy()
            while copy_graph != {}:
                path = []
                # Find semi-balanced start vertex
                for k,v in balance_dict.items():
                    if v[0]-v[1] == -1:
                        if k in copy_graph:
                            start_vertex = k
                        else:
                            # If start vertex not in graph anymore,
                            # Find a random start vertex
                            start_vertex = random.choice(list(get_vertices(copy_graph)))
                # Find initial edge
                vertex = start_vertex
                path.append(vertex)
                choices = copy_graph[vertex]
                original_vertex = vertex
                vertex = random.choice(choices)
                choices.remove(vertex)
                if choices == []:
                    del copy_graph[original_vertex]
                else:
                    copy_graph[original_vertex] = choices
                # Iterate untill current vertex == start vertex
                while vertex in copy_graph:
                    path.append(vertex)
                    original_vertex = vertex
                    choices = copy_graph[vertex]
                    vertex = random.choice(choices)
                    choices.remove(vertex)
                    if choices == []:
                        del copy_graph[original_vertex]
                    else:
                        copy_graph[original_vertex] = choices
                path.append(vertex)
                paths.append(path)
                vertices = list(get_vertices(copy_graph))
                for vertex in vertices:
                    if copy_graph[vertex] == []:
                        del copy_graph[vertex]
            print("The following paths were found:")
            for path in paths:
                print(path)
            print("\n")
            merge_paths(paths)
            # Merge paths
            total_path = paths[0]
            for path in paths[1:]:
                path = path[:-1]
                insert_pos = ''
                while insert_pos == '':
                    for i, vertex in enumerate(path):
                        if vertex in paths[0]:
                            path = path[i:]+path[:i]
                            insert_pos = paths[0].index(vertex)
                            break
                total_path[insert_pos:insert_pos] = path
            print("Merged they give this path:")
            print(total_path)
    else:
        print("Found cycle, not calculating eulerian trail.")

def merge_paths(paths):
    # Merge paths
    total_path = paths[0]
    for path in paths[1:]:
        path = path[:-1]
        insert_pos = ''
        while insert_pos == '':
            for i, vertex in enumerate(path):
                if vertex in paths[0]:
                    path = path[i:]+path[:i]
                    insert_pos = paths[0].index(vertex)
                    break
        total_path[insert_pos:insert_pos] = path
    print("Merged they give this path:")
    print(total_path)


def eulerian_wrapper(graph):
    # Wrapper function for assignment
    if assert_dict(graph):  # Check if graph is right format
        balance_dict = get_balance(graph)
        print(balance_dict)

        # If eulerian(balanced vertices)...
        if is_eulerian(balance_dict):
            print("Graph is balanced.")
            # ...find the cycle.
            find_eulerian_cycle(graph, balance_dict)

        # If vertices not balanced...
        else: 
            print("Unbalanced graph.")
            # ...check for eulerian path presence...
            if has_eulerian_path(graph, balance_dict):
                print("Found an eulerian path.")
                # ...and print this path.
                find_eulerian_path(graph, balance_dict)
        print("\n\n")
    else:
        print("Graph not a dictionary, exiting...")


if __name__ == "__main__":

    # GRAPH FROM FIG 8.22
    # Connected
    graph_822 = {'A':['B'],'B':['C'],'I':['H'],'H':['F'],'F':['G','E'],\
        'C':['I','J'],'G':['A'],'E':['J'],'J':['F','D'],'D':['C']}

    # Eulerian path, non-cycle
    graph_823 = {'A':['B'],'B':['C'],'I':['H'],'H':['F'],'F':['G','E'],\
        'C':['I','J'],'E':['J'],'J':['F','D'],'D':['C'], 'G':['K']}

    # Disconnected
    #graph_822 = {'A':['B'],'B':['C'],'I':['H'],'H':['F'],'F':['G','E'],\
        #'C':['I','J'],'E':['J'],'J':['F','D'],'D':['C'], 'G':['K'], 'Q':['Z']}

    # A SLIGHTLY BIGGER GRAPH, NEEDED FOR Q8
    bigger_graph = {1:[2], 2:[3], 3:[7,1],\
        4:[5,10],5:[6],6:[7],7:[8,9,11],\
        8:[4],9:[3,4],\
        10:[9],11:[12],12:[7]}
    # SPECTRUM FROM FIG 8.20
    s = ['ATG','TGG','TGC','GTG','GGC','GCA','GCG','CGT']

    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file) 
    eulerian_wrapper(graph_822)
    eulerian_wrapper(bigger_graph)
    eulerian_wrapper(graph_823)