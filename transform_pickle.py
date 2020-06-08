import os



if __name__=="__main__":

    ##nxpickled import
    dir = os.getcwd() + "/data/modules/CaRNAval/"

    try:
        import sys
        sys.path.append(os.path.abspath(dir))
        import RIN

    except:
        print("File not found : " + dir + "RIN.py")

    else:
        filename = "CaRNAval_1_as_dictionnary.nxpickled"

        try:
            import networkx
            import pickle

            objects = []

            with (open(dir+filename, "rb")) as openfile:
                while True:
                    try:
                        objects.append(pickle.load(openfile))
                    except EOFError:
                        break

            print("Dataset loaded")

        except OSError:
            print("File not found : " + dir + filename)

        else:


            ##Creation of a file for each RIN
            try:
                os.mkdir(dir + "Subfiles")
            except OSError:
                print ("Creation of the directory %s failed" % (dir + "Subfiles") + " : maybe it already exists ?")
            else:
                print ("Successfully created the directory %s " % (dir + "Subfiles"))

            header = "pos;k;seq\n"

            for i in range(1,338):
                motif = objects[0][i].graph
                f = open( dir + "Subfiles/" + str(i-1) + ".txt"  ,  "w+" )
                f.write(header)

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

                #print(components)

                num_nt = -1
                for a in components:
                    seq = ""
                    data_comp = str(num_nt+1)
                    for b in a:
                        num_nt += 1

                        #sometimes in the nxpicled file, a node has the attribute "realnt", and sometimes "real_nt", but it's the same thing
                        try:
                            seq += motif.nodes[b]["realnt"]

                        except:
                            seq += motif.nodes[b]["real_nt"]

                    data_comp += "," + str(num_nt) + ";" + str(len(a)) + ";" + seq + "\n"

                    f.write(data_comp)

                f.close()

                print(str(i-1) + ".txt created")

            print("Successfully parsed "+dir+filename)

