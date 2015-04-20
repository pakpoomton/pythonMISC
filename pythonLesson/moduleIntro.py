def sum_square(a,b):
    c = a*a + b*b
    return c

def sum_cube(a,b):
    c = a*a*a + b*b*b
    return c

def linear_seq(N):
    seq = [1]
    for x in range(0,N):
        ind = len(seq)
        last_element = seq[ind-1]
        seq = seq + [last_element + 1]
    return seq
