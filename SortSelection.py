def test():

    k       = 1 # k = 1 pra cada iteracao
    a       = np.array([[9,9,9,9,9,9],[1,2,3,4,5,6],[0,10,1,0,0,1],[2,3,4,5,6,7],[0,1,2,9,6,8],[1,2,3,4,5,6],[3,9,7,7,7,9]],np.float64)
    S       = np.full((7,7),-1,dtype=np.int32) # nao precisa setar
    d       = np.zeros(7,dtype=np.int32)       # seta zero pra cada iteracao
    sfront  = np.zeros(7,dtype=np.int32)       # nao precisa setar
    front   = np.zeros(7,dtype=np.int32)       # nao precisa setar
    rank    = np.zeros(7,dtype=np.int32)       # nao precisa setar
    result  = np.zeros((7,7),dtype=np.int32)
    density = np.zeros(7,dtype=np.float64)

    # Set true or false whether particle i dominates or not j
    for i in range(7):
        result[i]    = np.all( a[i] >= a, axis = 1)*np.any(a[i] > a, axis = 1)

    for i in range(7):

        # Find the particles that i dominates
        temp1 = np.where(result[i] == True)[0]
        # Put this particles in the set(i)
        S[i,:temp1.size] = temp1
        # End of set
        S[i][temp1.size] = -1
        # Count how many particles that dominate i
        d[i] = np.sum(result[:,i] == True)

    # Find all particles that are not dominated by anyone
    temp2 = np.where(d == 0)[0]
    # Put this particles in the pareto front
    front[:temp2.size] = temp2
    # Rank the pareto front's particles to 1
    rank[temp2] = 1

    front[temp2.size] = -1
    count    = temp2.size

    while count:

        i = count = 0

        while front[i] != -1:

            x = front[i]
            j = 0

            while S[x][j] != -1:

                y     = S[x][j]
                d[y] -= 1
                j    += 1

                if d[y] == 0:
                    rank[y]         = k + 1
                    sfront[count]   = y
                    count           += 1
            i += 1

        k               += 1
        sfront[count]   = -1
        front           = np.copy(sfront)

        #rank[rank == 0] = 10


    for i in range(2):
        x = a[:,i].argsort()
        b1 = (a[x[2:],i] - a[x[0:5],i] + 1)
        b2 = ((a[x[6],i] - a[x[0],i] + 1))

        b1*=b1
        b2*=b2

        density[x[1:6]] += b1 / b2
        density[x[0]] += 1.0
        density[x[6]] += 1.0

    print result
    print 'rank', rank
    print density
[ 1.49  0.08  2.    0.73  1.04  0.08  0.68]
