
import sys
import json
import random

input_key = 'Input'
output_key = 'OUT'
filter_key = 'Filter%d'

class NodeDAG:
    def __init__(self, f, inputs):
        self.f = f
        self.inputs = list(inputs)
    def __repr__(self):
        return 'NodeDAG(%r, %r)' % (self.f, self.inputs)
        
def checked_dag(dag):
    for i in range(len(dag)):
        if len(dag[i].inputs) != 1:
            assert dag[i].f[0] == 'FilterAdd', (dag, i, dag[i])
        for input in dag[i].inputs:
            if input >= i:
                print 'input >= i'
                print (i, dag[i].f, dag[i].inputs)
                print dag
                assert False
    return dag

def dag_to_json(dag):
    def sym(idx):
        if idx == -1:
            return input_key
        elif idx == len(dag)-1:
            return output_key
        else:
            return filter_key%idx
            
    ans = {input_key: ['ImageParam', 0]}
    checked_dag(dag)
    for i in range(len(dag)):
        ans[sym(i)] = [dag[i].f[0]] + [sym(x) for x in dag[i].inputs] + dag[i].f[1:]
    return ans
        
def to_pareto(dag_list):
    return [{'program': dag_to_json(dag), 'time': 0.0, 'error': 0.0, 'name': 'program%04d'%i} for (i, dag) in enumerate(dag_list)]

def write_pareto(pareto, filename):
    with open(filename, 'wt') as f:
        f.write(json.dumps(to_pareto(pareto), sort_keys=True, indent=4, separators=(',', ': ')))
    print 'Wrote %s' % filename

def node_has_recursion(node):
    return any(isinstance(y, NodeDAG) for y in node.f)

def node_shift_inputs(node, delta):
    return NodeDAG(node.f, [x+delta for x in node.inputs])
    
def flatten(dag):
    ans = []
    old_to_new = {-1: -1}
    for (i, node) in enumerate(dag):
        if node_has_recursion(node):
            len0 = len(ans)
            for sub in node.f:
                ans.append(node_shift_inputs(sub, len0))
        else:
            inputs = [old_to_new[x] for x in node.inputs]
            ans.append(NodeDAG(node.f, inputs))
        old_to_new[i] = len(ans)-1
    #print 'old_to_new:', old_to_new
    return ans
    
def main():
    args = sys.argv[1:]
    if len(args) < 1:
        print >> sys.stderr, 'topologies.py all'
        print >> sys.stderr, 'topologies.py initial imagew out_pareto.json filter_w'
        sys.exit(1)
    mode = args[0]
    
    filter_w = 9
    def append_iir(dag, input, xdir, ydir):
        if input is None:
            input = len(dag)-1
        K = [1.0 for i in range(filter_w)]
        F = [1.0 for i in range(filter_w)]
        return checked_dag(dag + [NodeDAG(['IIRFilter', K, F, xdir, ydir], [input])])

    def append_fir(dag, input=None, shape=None):
        if shape is None:
            shape = (fir_filter_w, fir_filter_w)
        if input is None:
            input = len(dag)-1
        K = [[random.random() for x in range(shape[1])] for y in range(shape[0])]
        return checked_dag(dag + [NodeDAG(['FIRFilter', K], [input])])
    
    def append_add(dag, inputL):
        return checked_dag(dag + [NodeDAG(['FilterAdd'], inputL)])

    def append_iir2_add(dag, dx=1, dy=0):
        input = len(dag)-1
        dag = append_iir(dag, input,  1*dx,  1*dy)
        dag = append_iir(dag, input, -1*dx, -1*dy)
        dag = append_add(dag, [input+1, input+2])
        return dag

    def append_fir2(dag, dx=1, dy=1, input=None):
        dag = append_fir(dag, len(dag)-1 if input is None else input, (1, fir_filter_w))
        dag = append_fir(dag, len(dag)-1,                             (fir_filter_w, 1))
        return dag

    def append_iir2(dag, dx=1, dy=1, input=None):
        dag = append_iir(dag, len(dag)-1 if input is None else input, dx,  0)
        dag = append_iir(dag, len(dag)-1,  0, dy)
        return dag

    def append_iir4(dag):
        dag = append_iir2_add(dag,  1,  0)
        dag = append_iir2_add(dag,  0,  1)
        return dag

    def multiband(dag_list):
        if len(dag_list) == 1:
            return dag_list[0]
        else:
            #print 'multiband: len(dag)=%d' % len(dag)
            dag1 = multiband(dag_list[1:])
            return flatten([NodeDAG(dag_list[0], [-1]), NodeDAG(['FilterDownsamplePrefilter'], [-1]), NodeDAG(dag1, [1]),
                           NodeDAG(['FilterUpsamplePostfilter'], [2]), NodeDAG(['FilterAdd'], [0, 3])])

    def append_svd(dag, n, iir):
        out = []
        for i in range(n):
            if iir:
                dag = append_iir4(dag)
            else:
                dag = append_fir2(dag, input=-1)
            out.append(len(dag)-1)
        if n > 1:
            dag = append_add(dag, out)
        return dag
        
    def append_iir3(dag, dx=1, dy=1):
        dag = append_iir(dag, len(dag)-1, dx,  0)
        dag = append_iir(dag, len(dag)-1,  0, dy)
        dag = append_iir(dag, len(dag)-1, dx, dy)
        return dag
            
    def append_iir4_add(dag):
        dag = append_iir2_add(dag, 1, 0)
        dag = append_iir2_add(dag, 0, 1)
        return dag
    
    if mode == 'all':
        fir_filter_w = 30*2+1 #31
        max_levels = 7
        
        write_pareto([append_iir4_add([])], 'topologies_iir4.json')
        write_pareto([append_iir2([])], 'topologies_iir2.json')
        for n in range(1,max_levels+1):
            write_pareto([append_svd([], n, True )], 'topologies_svd_iir%d.json'%n)
            write_pareto([append_svd([], n, False)], 'topologies_svd_fir%d.json'%n)
            
        write_pareto([append_iir3([])], 'topologies_iir3.json')
        write_pareto([append_iir3(append_iir3([]))], 'topologies_iir3_twice.json')
        write_pareto([append_iir2(append_fir([]))], 'topologies_iir2_fir.json')
        write_pareto([append_iir2(append_iir2([]))], 'topologies_iir2_twice.json')
    
        fir_filter_w = 9 #15
        for bands in range(1, 10):
            for n in range(1,max_levels+1):
                write_pareto([multiband([append_svd([], n, False)]*bands)], 'topologies_bands%d_fir%d.json'%(bands, n))
                fir_level = append_svd([], n, False)
                iir_level = append_svd([], n, True )
                write_pareto([multiband([fir_level]*(bands-1)+[iir_level])], 'topologies_bands%d_iir%d.json'%(bands, n))
    #            write_pareto([multiband([iir_level]+[fir_level]*(bands-1))], 'topologies_bands%d_iir%d.json'%(bands, n))
    elif mode == 'initial':
        fir_filter_w = int(args[3]) #9 #15
        imagew = int(args[1])
        max_bands = 1
        w_eff = fir_filter_w
        while w_eff < imagew:
            w_eff *= 2
            max_bands += 1
        
        ans = []
        for bands in range(1, max_bands+1):
            ans.append(multiband([append_svd([], 1, False)]*bands))
        out_pareto = args[2]
        print 'Writing Pareto to %s (max_bands=%d, filter_w=%d)' % (out_pareto, max_bands, fir_filter_w)
        write_pareto(ans, out_pareto)
    else:
        raise ValueError('bad mode %s'%mode)
    
if __name__ == '__main__':
    main()
    
