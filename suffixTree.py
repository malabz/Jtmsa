'''
build suffix tree with ukk algo in O(n) time

author: Juntao Chen
date: 12.29.2021
'''

end = -1

class Node(object):
    def __init__(self, leaf:bool):
        # the start & end
        self.start = None
        self.end = None
        self.leaf = leaf
        # self.length = 0
        self.children = {}
        self.suffix_link = None
    
    def __getattribute__(self, name:str):
        if name == 'end':
            if self.leaf:
                return end
        return super(Node, self).__getattribute__(name)

class suffixTree(object):
    def __init__(self, strings:str):
        self.T = strings + '#'
        self.actNode = None
        self.actEdge = -1
        self.actLength = 0
        self.remainder = 0
        self.insideNode = None
        # the length of strings
        self.len = -1
        self.root = None
        self.build_tree()

    def _edge_length(self, node:Node):
        return node.end - node.start + 1

    def _walk_down(self, node:Node):
        length = self._edge_length(node)
        if self.actLength >= length:
            self.actEdge += length
            self.actLength -= length
            self.actNode = node
            return True
        return False

    def _gen_node(self, start:int, end=None, leaf=False):
        node = Node(leaf)
        node.suffix_link = self.root
        node.start = start
        # end is None which means it is a leaf node.
        node.end = end
        return node
    
    def _gen_trie(self, pos:int):
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
                self.actNode.children[self.T[self.actEdge]] = self._gen_node(pos, leaf=True)
                if self.insideNode:
                    self.insideNode.suffix_link = self.actNode
                    self.insideNode = None

            else:
                nextNode = self.actNode.children.get(self.T[self.actEdge])
                '''当前节点是否存在跨点'''
                if self._walk_down(nextNode):
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
                split_node = self._gen_node(nextNode.start, splitEnd)
                # 活跃点的孩子重新设置
                self.actNode.children[self.T[self.actEdge]] = split_node
                # 将两个孩子节点添加至内部节点上
                split_node.children[self.T[pos]] = self._gen_node(pos, leaf=True)
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
    
    def build_tree(self):
        self.len = len(self.T)
        rootEnd = -1
        self.root = self._gen_node(-1, rootEnd)
        self.actNode = self.root
        for i in range(self.len): self._gen_trie(i)