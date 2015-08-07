
import copy

FUNCNAME = 0
INPUTS = 1

def simplify(p):
    ans = copy.deepcopy(p)
    seen = set()
    def convert_recursive(image):
        if id(image) in seen:
            return
        seen.add(id(image))
        if image[FUNCNAME] == '+':
            while True:
                changed = False
                i = 0
                while i < len(image[INPUTS]):
                    input = image[INPUTS][i]
                    if input[FUNCNAME] == '+':
                        image[INPUTS][i:i+1] = input[INPUTS]
                        changed = True
                    i += 1
                if not changed:
                    break
        image[INPUTS] = sorted(image[INPUTS], key=lambda x: x[0])  # Sort comparison key
#        image[INPUTS] = sorted(image[INPUTS])                     # This works equally well
        
    convert_recursive(ans)
    return ans
    
def program1():
    I1 = ['FilterFIR(p1)', []]
    I2 = ['+', [I1, I1]]
    return I2

def program2():
    I1 = ['FilterFIR(p1)', []]
    I2 = ['FilterFIR(p2)', []]
    I3 = ['+', [I1, I2]]
    return I3

def program3():
    I1 = ['FilterFIR(p1)', []]
    I2 = ['FilterFIR(p2)', []]
    I3 = ['+', [I2, I1]]
    return I3

def program4():
    I1 = ['FilterFIR(p1)', []]
    I2 = ['FilterFIR(p2)', []]
    I3 = ['+', [I1, I1]]
    I4 = ['+', [I2, I1]]
    return I4

def program5():
    I1 = ['FilterFIR(p1)', []]
    I2 = ['FilterFIR(p2)', []]
    I3 = ['+', [I2, I1]]
    I4 = ['+', [I3, I1]]
    return I4

def program6():
    I1 = ['FilterFIR(p1)', []]
    I2 = ['FilterFIR(p2)', []]
    I3 = ['+', [I1, I2]]
    I4 = ['+', [I1, I3]]
    return I4

def program7():
    I1 = ['FilterFIR(p1)', []]
    I2 = ['FilterFIR(p2)', []]
    I3 = ['+', [I1, I1]]
    I4 = ['+', [I3, I1]]
    I5 = ['+', [I4, I2]]
    return I5

def program8():
    I1 = ['FilterFIR(p1)', []]
    I2 = ['FilterFIR(p2)', []]
    I3 = ['+', [I2, I1]]
    I4 = ['+', [I3, I1]]
    I5 = ['+', [I4, I1]]
    return I5

def test():
    sp1 = simplify(program1())
    sp2 = simplify(program2())
    sp3 = simplify(program3())
    sp4 = simplify(program4())
    sp5 = simplify(program5())
    sp6 = simplify(program6())
    sp7 = simplify(program7())
    sp8 = simplify(program8())
    
    assert sp1 != sp2
    assert sp3 == sp2
    assert sp4 == sp3
    assert sp5 != sp4
    assert sp5 != sp3
    assert sp5 != sp2
    assert sp5 != sp1
    assert sp6 == sp5
    assert sp7 not in [sp1, sp2, sp3, sp4, sp5, sp6]
    assert sp8 == sp7
    print 'simplify:   OK'

def main():
    test() #print simplify(program1())

if __name__ == '__main__':
    main()
    
