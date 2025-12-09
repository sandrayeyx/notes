#segmath运行
ciphertext=[62,117,16,
            108,57,10,
            112,73,44,
            48,4,94,
            74,0,105,
            98,3,12,
            21,56,100,
            31,8,73,
            120,107,10,
            70,25,102,
            122,20,80,
            84,4,3,
            62,82,126]
R = Zmod(127)
vcode = {}
for i in range(len(ciphertext)//3):
    vcode[i]=vector([R(ciphertext[3*i]),R(ciphertext[3*i+1]),R(ciphertext[3*i+2])])
print(vcode)
def solve1():
    for guess1 in range(0x20,127):
        for guess2 in range(0x20,127):
            mc = Matrix(R,3,3,[ord('v'),ord('m'),ord('c'),
                              ord('{'),guess1,guess2,
                              0x20,0x20,0x20])
            me = Matrix(R,3,3,[ciphertext[0],ciphertext[1],ciphertext[2],
                               ciphertext[3],ciphertext[4],ciphertext[5],
                               ciphertext[-3],ciphertext[-2],ciphertext[-1]])
            try:
                key = mc.solve_right(me)
                cipher = ''
                key = key.inverse()
                for i in range(len(ciphertext)//3):
                    v = vcode[i]*key
                    for j in range(3):
                        assert (v[j]==0x20) or chr(v[j]).isalnum() or chr(v[j]) in '{}_'
                    cipher = cipher+chr(v[0])+chr(v[1])+chr(v[2])
                print(cipher)
                return 1
            except:
                pass
    return 0
def solve2():
    for guess1 in range(0x20,127):
        if (guess1 == 0x20) or chr(guess1).isalnum() or chr(guess1) in '{}_':
            for guess2 in range(0x20,127):
                if (guess2 == 0x20) or chr(guess2).isalnum() or chr(guess2) in '{}_':
                    for guess3 in range(0x20,127):
                        if (guess3 == 0x20) or chr(guess3).isalnum() or chr(guess3) in '{}_':
                            mc = Matrix(R,3,3,[ord('v'),ord('m'),ord('c'),
                                               ord('{'),guess1,guess2,
                                               guess3,ord('}'),0x20])
                            me=Matrix(R,3,3,[ciphertext[0],ciphertext[1],ciphertext[2],
                                             ciphertext[3],ciphertext[4],ciphertext[5],
                                             ciphertext[-3],ciphertext[-2],ciphertext[-1]])
                        try:
                            key = mc.solve_right(me)
                            cipher = ''
                            key=key.inverse()
                            for i in range(len(ciphertext)//3):
                                v=vcode[i]*key
                            for j in range(3):
                                assert (v[j]==0x20) or chr(v[j]).isalnum() or chr(v[j]) in '{}_'
                            cipher = cipher+chr(v[0])+chr(v[1])+chr(v[2])
                            print(cipher)
                            return 1
                        except:
                            pass
    return 0
def solve3():
    for guess1 in range(0x20,127):
        for guess2 in range(0x20,127):
            mc = Matrix(R,3,3,[ord('v'),ord('m'),ord('c'),
                              ord('{'),guess1,guess2,
                              ord('}'),0x20,0x20])
            me=Matrix(R,3,3,[ciphertext[0],ciphertext[1],ciphertext[2],
                            ciphertext[3],ciphertext[4],ciphertext[5],
                            ciphertext[-3],ciphertext[-2],ciphertext[-1]])
            try:
                key = mc.solve_right(me)
                cipher = ''
                key=key.inverse()
                for i in range(len(ciphertext)//3):
                    v=vcode[i]*key
                    for j in range(3):
                        assert ((v[j]==0x20) or chr(v[j]).isalnum() or chr(v[j]) in '{}_')
                    cipher = cipher+chr(v[0])+chr(v[1])+chr(v[2])
                print(cipher)
                return 1
            except:
                pass
    return 0
a=solve1()
print(a)
a=solve3()
print(a)
a=solve2()
print(a)