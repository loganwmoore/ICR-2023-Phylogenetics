from ete3 import PhyloTree, TreeStyle


t = PhyloTree("(((A:2,A:2),B:3),B:4);")



#def sumstat1(tree):
    
#def sumstat2(tree):




#def sumstat3(tree):



def sumstat5(tree):
    numcount = 0
    firstA = tree.search_nodes(name="A")[0]
    total = 0
    name = tree.get_leaves_by_name("A")
    for n in name:
        numcount+=1
        total+= firstA.get_distance(n)
    numcount-=1
    return total/numcount
#def sumstat6(tree):




#def sumstat7(tree):
    
print(sumstat5(t))
t.show()
