#!/usr/bin/env python3

import sys
# import networkx as nx
# import matplotlib.pyplot as plt
# import pydot
# from networkx.drawing.nx_pydot import graphviz_layout

def preText(s):
    output = ""
    for x in s:
        if x == '_':
            output += '\\_'
        else:
            output += x
    return output

def processNode(graph_stack, labels, level, graph):
    todo = []
    while(graph_stack[-1] > level):
        graph_stack.pop()
        todo.append(labels.pop())
    # print(labels[-1], "->")
    for x in todo:
        # print(x)
        graph.add_edge(labels[-1], x)
    # print("")

if __name__ == "__main__":
    file = open(sys.argv[1], 'r')

    total_times = {}
    graph_stack = [-1]
    # labels = [None]
    last_function = ""
    last_level = -1
    c = 0
    # function_calls = nx.DiGraph()

    while True:
        line = file.readline()

        if not line:
            break

        line = line.split()

        if (len(line) == 3):
            if line[1] not in total_times:
                total_times[line[1]] = 0
            if line[0] == "Start:":
                level = int(line[2])
                # if(graph_stack[-1] > level):
                    # processNode(graph_stack, labels, level, function_calls)
                graph_stack.append(level)
                c += 1
                #labels.append(line[1]+'_'+str(c))
                ## labels.append(line[1])
            if line[0] == "End:":
                last_level = int(line[2])
                last_function = line[1]

        if (len(line) == 2):
            total_times[last_function] += float(line[1])

            if last_level == 0:
                sorted_total_times = dict(sorted(total_times.items(), key=lambda item: -item[1]))
                for func_name, time in sorted_total_times.items():
                    print("\\texttt{", preText(func_name), "} &", round(time, 3),"\\\\")
                print()
                # processNode(graph_stack, labels, 0, function_calls)
                # pos = graphviz_layout(function_calls, prog="twopi", root="find_certificate_0")
                # nx.draw(function_calls, pos, with_labels=True, arrows=True)
                # plt.show()

                # total_times = {}
                # graph_stack = [-1]
                ## labels = [None]
                # last_function = ""
                # last_level = -1
                # c = 0
                # function_calls = nx.DiGraph()

                graph_stack.pop()
                # labels.pop()
