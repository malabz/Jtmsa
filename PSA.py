from PSA_LCS import PSA_LCS
from PSA_Kband import PSA_AGP_Kband
from PSA_Suffixtree import PSA_Suffixtree

class PSA(object):
    def __init__(self, A, B, mode="suffix"):
        self._state_loc = 0
        if len(A) >= len(B):
            self._state_loc = 1
            A, B = B, A
        self.A = A
        self.B = B
        self.mode = mode
    
    def _build_dict(self):
        dict_A, dict_B, dict_ = {}, {}, {}
        for i in self.A:
            if dict_A.get(i) is None:
                dict_A[i] = 1
            else:
                dict_A[i] += 1
        for i in self.B:
            if dict_B.get(i) is None:
                dict_B[i] = 1
            else:
                dict_B[i] += 1
        len_A = len(dict_A)
        len_B = len(dict_B)
        dict_ = dict_A.copy()
        dict_.update(dict_B)
        len_ = len(dict_)

        self.union = False if len_ >= len_A + len_B else True
        
        return dict_A, dict_B, dict_
    
    def align(self):
        tree_A = PSA_Suffixtree(self.A)
        return tree_A.align(self.B)

