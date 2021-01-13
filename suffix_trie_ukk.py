# -*- coding: utf-8 -*-
# build suffix trie with ukk algo in O(n) time

import networkx as nx
import matplotlib.pyplot as plt
from pprint import pprint

end = -1

class Node(object):
    def __init__(self, leaf):
        # the start & end
        self.start = None
        self.end = None
        self.leaf = leaf
        self.children = {}
        self.suffix_link = None
    
    def __getattribute__(self, name: str):
        if name == 'end':
            if self.leaf:
                return end
        return super(Node, self).__getattribute__(name)

class Suffix_Trie(object):
    def __init__(self, strings):
        self.T = strings
        self.actNode = None
        self.actEdge = -1
        self.actLength = 0
        self.remainder = 0
        self.insideNode = None
        # the length of strings
        self.len = -1
        self.root = None

    def edge_length(self, node):
        return node.end - node.start + 1

    def walk_down(self, node):
        length = self.edge_length(node)
        if self.actLength >= length:
            self.actEdge += length
            self.actLength -= length
            self.actNode = node
            return True
        return False

    def gen_node(self, start, end=None, leaf=False):
        node = Node(leaf)
        node.suffix_link = self.root
        node.start = start
        # end is None which means it is a leaf node.
        node.end = end
        return node
    
    def gen_trie(self, pos):
        global end
        # init
        end = pos
        self.remainder += 1
        self.insideNode = None

        while self.remainder > 0:
            if self.actLength == 0:
                self.actEdge = pos
            
            if self.actNode.children.get(self.T[self.actEdge]) is None:
                # creat the leaf node
                self.actNode.children[self.T[self.actEdge]] = self.gen_node(pos, leaf=True)
                if self.insideNode:
                    self.insideNode.suffix_link = self.actNode
                    self.insideNode = None

            else:
                nextNode = self.actNode.children.get(self.T[self.actEdge])
                '''当前节点是否存在跨点'''
                if self.walk_down(nextNode):
                    continue
                '''当前字符在边上'''
                if self.T[nextNode.start + self.actLength] == self.T[pos]:
                    if self.insideNode and (self.actNode != self.root):
                        self.insideNode.suffix_link = self.actNode
                        self.insideNode = None

                    self.actLength += 1
                    break
                '''当前字符不在边上'''
                splitEnd = nextNode.start + self.actLength - 1
                # 新的内部节点
                split_node = self.gen_node(nextNode.start, splitEnd)
                # 活跃点的孩子重新设置
                self.actNode.children[self.T[self.actEdge]] = split_node
                # 将两个孩子节点添加至内部节点上
                split_node.children[self.T[pos]] = self.gen_node(pos, leaf=True)
                # 重新设置节点的范围
                nextNode.start += self.actLength
                split_node.children[self.T[nextNode.start]] = nextNode

                if self.insideNode:
                    self.insideNode.suffix_link = split_node
                self.insideNode = split_node

            self.remainder -= 1
            if (self.actNode == self.root) and (self.remainder > 0):
                self.actLength -= 1
                self.actEdge = pos - self.remainder + 1
            elif self.actNode != self.root:
                self.actNode = self.actNode.suffix_link
    
    def suffix_trie(self):
        self.len = len(self.T)
        rootEnd = -1
        self.root = self.gen_node(-1, rootEnd)
        self.actNode = self.root
        for i in range(self.len):
            self.gen_trie(i)
        
    def creat_graph(self, node, deepth=1, edges=[]):
        if node.start != -1:
            strs = self.T[node.start:(node.end+1)]
        else:
            strs = 'Root'
        if node.children:
            sons = node.children.keys()
            for k, son in enumerate(sons):
                strs_son = T[node.children[son].start:(node.children[son].end+1)]
                edges.append([strs, strs_son, deepth, k])
                self.creat_graph(node.children[son], deepth+1, edges)
        
        return edges
    
    def draw(self):
        edges = self.creat_graph(self.root)
        G = nx.DiGraph()
        # pprint(edges)

        i, j = 0, 1
        pos = {}
        pos[i] = (0,0)
        record = {0:{}}
        for k in range(len(edges)):
            print(edges[k][0]+'-->'+edges[k][1])
            record[edges[k][2]] = {}
            if k >= 1 :
                if edges[k][0] == 'Root':
                    i = 0
                else:
                    i = record[edges[k][2]-1][edges[k][0]]
                j += 1
            record[edges[k][2]][edges[k][1]] = j
            pos[j] = (edges[k][3]+2*k, -edges[k][2])
            G.add_edge(i,j, name=edges[k][1])

        nx.draw(G, pos)
        edge_labels = nx.get_edge_attributes(G, 'name')
        nx.draw_networkx_edge_labels(G,pos, edge_labels=edge_labels)
        plt.show()
        
if __name__ == "__main__":
    T = "GGGACGGCGGGACGCAGGGA#"
    tree = Suffix_Trie(T)
    tree.suffix_trie()
    tree.draw()





        

