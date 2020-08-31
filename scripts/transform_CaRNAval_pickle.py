#!/usr/bin/python3

# This script's purpose is to extract information about the CaRNAval
# RINS from a Python pickle object containing RINs from their RIN.py class.
# We do this because the official JSON file is hard to understand, and Antoine Soulé
# recommended the pickle.

import networkx, os, pickle, sys

if __name__=="__main__":

    
    rin_DIR = os.getcwd() + "/../data/modules/RIN/"
    filename = "CaRNAval_1_as_dictionnary.nxpickled"

    # Check that we can find CaRNAval RINs, and load the dataset
    try:
        sys.path.append(os.path.abspath(rin_DIR))
        import RIN
    except:
        print("File not found:" + rin_DIR + "RIN.py")
        exit(1)

    try:
        objects = []
        with (open(rin_DIR+filename, "rb")) as openfile:
            while True:
                try:
                    objects.append(pickle.load(openfile))
                except EOFError:
                    break
        print("Dataset loaded")
    except OSError:
        print("File not found : " + rin_DIR + filename)
        exit(1)

    # Creation of a directory to extract RINs from the pickle file to individual files
    try:
        os.makedirs(rin_DIR + "Subfiles", exist_ok=True)
    except OSError:
        print("Creation of the directory %s failed" % (rin_DIR + "Subfiles"))
        exit(1)

    # Loop on every CaRNAval module and extract it from the Python object to flat text file
    n_modules = len(objects[0]) # ? to
    for i in range(1,1+n_modules):
        motif = objects[0][i].graph
        f = open(rin_DIR + "Subfiles/" + str(i-1) + ".txt", "w+")
        f.write("ntA,ntB,long_range;...\n")

        components = []
        comp = []
        nodes = list(motif)
        nodes.sort()
        for node in nodes:
            if comp == []:
                comp.append(node)
            else:
                if comp[-1] + 1 != node : #not the same component
                    components.append(comp)
                    comp = []
                    comp.append(node)
                else :
                    comp.append(node)
        components.append(comp)

        #print(nodes)

        basepairs = ""
        edges = list(motif.edges())
        for a in edges:
            if motif.edges[a]['label'] == 'CWW' :
                ntA = nodes.index(a[0])
                ntB = nodes.index(a[1])

                if ntA <= ntB :
                    basepairs += str(ntA) + "," + str(ntB) + "," + str(motif.edges[a]['long_range']) + ";"

        f.write(basepairs + "\n")
        f.write("pos;k;seq\n")

        num_nt = -1
        for a in components:
            seq = ""
            data_comp = str(num_nt+1)
            for b in a:
                num_nt += 1

                # sometimes in the nxpicled file, a node has the attribute "realnt", 
                # and sometimes "real_nt", but it's the same thing
                try:
                    seq += motif.nodes[b]["realnt"]
                except:
                    seq += motif.nodes[b]["real_nt"]
            data_comp += "," + str(num_nt) + ";" + str(len(a)) + ";" + seq + "\n"
            f.write(data_comp)

        f.close()
        # print(str(i-1) + ".txt created")

    print("Successfully parsed "+filename, ", now individual RINs are saved in Subfiles/ folder.", sep='')

