def qsort(a,lo,hi):
    if lo<hi:
        p=partition(a,lo,hi)
        qsort(a,lo,p-1)
        qsort(a,p+1,hi)
    return a


def partition(a,lo,hi):
    pivot=a[hi]
    i=lo
    for j in range(lo,hi):
        if a[j]<pivot:
            a[j],a[i]=swap(a[j],a[i])
            i +=1 
    a[i],a[hi]=swap(a[i],a[hi])
    return i

def swap(a,b):
    return b,a

a=[6,4,8,1]
print(a)
print(qsort(a,0,3))


